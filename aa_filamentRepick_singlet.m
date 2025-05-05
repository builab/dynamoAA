%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to apply alignment parameters to repick filament with torsion model
% Should have same parameters as imodModel2Filament
% Now successfully do initial angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Before Running Script %%%%%%%%%%
%% Activate Dynamo
%run /london/data0/software/dynamo/dynamo_activate.m
run /data2/apps/dynamo/1.1.546/dynamo_activate.m 

% Change path to the correct directory
prjPath = '/mnt/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/Singlet_STA/';


%%%%%%% Variables subject to change %%%%%%%%%%%
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
modelDir = sprintf('%smodels_repick', prjPath);
origParticleDir = sprintf('%sparticles', prjPath);
particleDir = sprintf('%sparticles_repick', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
pixelSize = 14.00; % Angstrom per pixel
periodicity = 83.2; % Measured from yor average
boxSize = 60;
mw = 12;
subunits_dphi = 0.0;  % For the tip CP 0.72, baseCP 0.5, doublet 0
subunits_dz = periodicity/pixelSize; % in pixel repeating unit dz = 8.4 nm = 168 Angstrom/pixelSize
filamentRepickListFile = sprintf('%sfilamentRepickList.csv', prjPath);
filamentListFile = sprintf('%sfilamentListPolarity.csv', prjPath);
tableAlnFileName = 'merged_particles_singlet_align.tbl'; % merge particles before particle alignment for robust but must be merged_particles_align to use doInitialAngle
avgLowpass = 30; % Angstrom
tomoSuffix = '_14.00Apx';
dTh = 30; % Distance Threshold in Angstrom
doExclude = 1; % Exclude particles too close
doOutlier = 1; % Exclude outlier using CC using MAD
doInitialAngle = 1; % Should turn on everything except pure microtubule

%%%%%%% Do not change anything under here %%%%%

%% loop through all tomograms
filamentList = readcell(filamentListFile, 'Delimiter', ',');
fileID = fopen(docFilePath); D = textscan(fileID,'%d %s'); fclose(fileID);
tomoID = D{1,1}'; % get tomogram ID
nTomo = length(D{1,2}); % get total number of tomograms

tAll = dread(tableAlnFileName);

filamentRepickList = {};
particleCounts = [];  % Array to store the number of particles per filament

%% Loop through tomograms
for idx = 1:nTomo
    tomo = D{1,2}{idx,1};
    [tomoPath,tomoName,ext] = fileparts(tomo);
    tomono = D{1,1}(idx);
    % Modify specific to name
    tomoName = strrep(tomoName, tomoSuffix, ''); % Remove the suffix part of the name
    tTomo = tAll(tAll(:,20) == tomono, :);
    if isempty(tTomo) == 1
        continue;
    end
      
    modelout =   [modelDir '/' tomoName '.omd'];
    contour = unique(tTomo(:, 23));
    
    m = {}; % Cell array contains all filament
    
    for i = 1:length(contour)
        tContour = tTomo(tTomo(:, 23) == contour(i), :);

        
	%%%% Eliminating outliner using mean absolute deviation      
        if doOutlier > 0
            cc = tContour(:, 10);
            x = median(cc);
            y = mad(cc);
            tContour = tContour(cc > x - 3*y, :);
            disp(['Contour ' num2str(contour(i)) ': Exclude ' num2str(sum(cc <= x - 3*y)) ' particles']);
        end
        %%%%

        %%%% Code to estimate the Phi angle of the middle section        
	% Check if there are at least 10 rows
	if size(tContour, 1) < 10
    		phi = median(tContour(:, 9));
	else
		% Determine the middle region indices
		middleStart = round((size(tContour, 1) - 10) / 2) + 1;  % Start index of middle 10 rows
		middleEnd = middleStart + 9;  % End index of middle 10 rows
		% Extract column 9 from the middle region and compute the median
		phi = median(tContour(middleStart:middleEnd, 9));
	end
	% This actually make less error to find
	[~, phi_idx] = min(abs(tContour(:, 9) - phi));
	phi_xyz = tContour(phi_idx, 24:26);
	fprintf('Median of Phi: %.4f\n', phi);
	%%%
    
        %%%% Exclude duplicate particles
        if doExclude > 0
            tContourEx = dpktbl.exclusionPerVolume(tContour, dTh/pixelSize);
            % Make sure to sort by particles number for not inverting angle
            tContour = sortrows(tContourEx, 1);
            disp(['Exclude ' num2str(size(tContour, 1) - size(tContourEx, 1)) ' particles due to proximity']);
        end
        %%%% End exclude
		
	%%%% Flip polarity if is shown in the filamentList
        filamentName = [tomoName '_' num2str(contour(i))];
        % Look up for polarity
        filIdx = find(strcmp(filamentList(:, 1), filamentName));
        polarity = filamentList{filIdx, 2};
        if polarity > 0
            tContour = flipud(tContour);
            disp([filamentName ' - flip polarity']);
	else
		disp([filamentName ' - no flip']);
	end
        %%%% End flipping
   
        if isempty(tContour) == 1
            continue;
        end

        points = tContour(:, 24:26) + tContour(:, 4:6);
	if size(points, 1) < 5
		disp('Less than 5 particles to fit');
		continue;
	end
      
        m{i} = dmodels.filamentWithTorsion();
        m{i}.subunits_dphi = subunits_dphi;
        m{i}.subunits_dz = subunits_dz;
        
        m{i}.name = [tomoName '_' num2str(contour(i))];
        % Import coordinate
        m{i}.points = points;
        % Create backbone
        m{i}.backboneUpdate();
        % Update crop point (can change dz)
        m{i}.updateCrop();
        % Link to catalog
        m{i}.linkCatalogue(c001Dir, 'i', idx);
        m{i}.saveInCatalogue();
        t = m{i}.grepTable();
        
        %%%% In case of empty
        if isempty(t) == 1
          	warning(['Skip: ' tomoName  '_' num2str(contour(i)) 'does not have any particles!']);
        	continue;
        end
        
        %%%% Set contourID
	fprintf('Assign contourID %d\n', i);
        t(:,23) = contour(i); % Additing contour number (filament)
        
        if doInitialAngle > 0
        	% Find the point in t nearest to phi_xyz
        	% Compute squared Euclidean distances to avoid unnecessary sqrt computation
        	distances = sum((t(:, 24:26) - phi_xyz).^2, 2);
		% Find the index of the closest point
		[~, minIndex] = min(distances);
            	% Set the phi angle of the minIndex same as phi
            	t(:, 9) = t(:, 9) - t(minIndex, 9) + phi; 
        end
        
        % Check point
        dwrite(t, [modelDir '/' tomoName '_' num2str(contour(i)) '.tbl']);
        targetFolder = [particleDir '/'  tomoName '_' num2str(contour(i))];
        
        % Cropping subtomogram out
        % 0.2b
        try
       		dtcrop(docFilePath, t, targetFolder, boxSize, 'mw', mw);
        	tCrop = dread([targetFolder '/crop.tbl']);
        	oa_all = daverage(targetFolder, 't', tCrop, 'fc', 1, 'mw', mw);
        	dwrite(dynamo_bandpass(oa_all.average, [1 round(pixelSize/avgLowpass*boxSize)]), [targetFolder '/average.em']);
        	% Average the middle region again
        	if size(tCrop, 1) > 15
            	midIndex = floor(size(tCrop, 1)/2);
            	tCrop = tCrop(midIndex - 3: midIndex + 4, :);
        	end
        	oa = daverage(targetFolder, 't', tCrop, 'fc', 1, 'mw', mw);
        	dwrite(dynamo_bandpass(oa.average, [1 round(pixelSize/avgLowpass*boxSize)]), [targetFolder '/template.em']);
        
        	% Plotting save & close. dtplot seems to error if only 1 particles
        	if size(tCrop, 1) > 1
            	dtplot(tCrop, 'pf', 'oriented_positions');
            	view(-230, 30); axis equal;
            	print([targetFolder '/repick_' tomoName '_' num2str(contour(i))] , '-dpng');
            	close all
        	end
       	catch
  			warning(['Skip: ' tomoName  '_' num2str(contour(i)) 'does not have enough particles!'])
  			continue;
        end
        % Add to particle count
        particleCounts(end + 1) = size(t, 1);
        
  		% If it is cropping out
  		filamentRepickList{end + 1, 1} = [tomoName  '_' num2str(contour(i))];   
    end
    % Write the DynamoModel
    dwrite(m, modelout);
end

%% Generate and Save Histogram of Particle Counts per Filament
if ~isempty(particleCounts)
    figure('Name', 'Particles per Filament', 'NumberTitle', 'off');
    histogram(particleCounts, 'BinMethod', 'auto');
    histogram(particleCounts, 10);  % 10 bins
    xlabel('Number of Particles per Filament');
    ylabel('Frequency');
    title('Histogram of Particle Counts per Filament');
    grid on;
    histogramFile = sprintf('%sparticleCountsHistogram.png', prjPath);
    saveas(gcf, histogramFile);
    disp(['Histogram saved to ' histogramFile]);
    
    %close(gcf);
else
    disp('No particle counts to plot.');
end

%% Write filament list out
writecell(filamentRepickList, filamentRepickListFile);


