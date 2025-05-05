%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to convert IMOD model to filament torsion model covering the missing wedge
% dynamoDMT v1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using GUI https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Filament_model
% Imod coordinate should be in text file, clicking along the filament (no direction needed)
% model2point -Contour imodModel.mod imodModel.txt
% Rationale: Generate particles with "8nm + rise" and rotate 360/13 according to 13PF


%%%%%%%% Before Running Script %%%%%%%%%%%%%%%
%%% Activate Dynamo
run /london/data0/software/dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/Singlet_STA/';

%%%%%%% Variables subject to change %%%%%%%%%%%

docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
modelDir = sprintf('%smodels', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
recSuffix = '_14.00Apx'; % The suffix path without .mrc
pixelSize = 14.00; % Angstrom per pixel
periodicity = 83.2; % Using 84.5 of doublet, 82.8 for CP tip, 86 for CP base
subunits_dphi = 0.0;  % For the tip CP 0.72, base CP 0.5, doublet 0
subunits_dz = periodicity/pixelSize;
filamentListFile = sprintf('%sfilamentList.csv', prjPath);
minPartNo = 2; % Minimum particles number per Filament

%%%%%%% Do not change anything under here %%%%%

% loop through all tomograms
fileID = fopen(docFilePath); D = textscan(fileID,'%d %s'); fclose(fileID);
tomoID = D{1,1}'; % get tomogram ID
nTomo = length(D{1,2}); % get total number of tomograms

filamentList = {};
particleCounts = [];  % Array to store the number of particles per filament

%% Loop through tomograms
for idx = 1:nTomo
    tomo = D{1,2}{idx,1};
    [tomoPath,tomoName,ext] = fileparts(tomo);
    % Modify specific to name
    tomoName = strrep(tomoName, recSuffix, ''); % Remove the rec part of the name from IMOD
    imodModel = [modelDir '/' tomoName '.txt'];
    modelout = strrep(imodModel, '.txt', '.omd');
    
    allpoints = load(imodModel);
    
    m = {}; % Cell array contains all filament
    contour = unique(allpoints(:, 1));
    % Loop through filaments
    for i = 1:length(contour)
        filamentid = contour(i);
        points = allpoints(allpoints(:, 1) == filamentid, 2:4);
        % Sort it so small Y first
        points = sortrows(points, 2);
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
        
        % Testing this block
        t = m{i}.grepTable();

        % 0.2b addition
        t(:,23) = contour(i);
        if (size(t, 1) < minPartNo)
        	disp(['Skip ' tomoName ' Contour ' num2str(contour(i)) ' with less than ' num2str(minPartNo) ' particles'])
        	continue
        end
        % Add to particle count
        particleCounts(end + 1) = size(t, 1);

        % Add the good to the list
        filamentList{end + 1, 1} = [tomoName '_' num2str(contour(i))];
		dwrite(t, [modelDir '/' tomoName '_' num2str(contour(i)) '.tbl']);
        
        % Optional for visualization of table
        dtplot(t, 'pf', 'oriented_positions');
        view(-230,30);axis equal;
        hold on;
    end
    print([modelDir '/' tomoName] , '-dpng');
    close all;
    
    % Write the DynamoModel
    dwrite(m, modelout)

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

%% Write out list file
writecell(filamentList, filamentListFile);
