%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align filament to a reference
% dynamoMT v0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The new list will contain 2 columns (Filament, polarity)

%%%%%%%% Before Running Script %%%%%%%%%%%%%%%

%%% Activate Dynamo
run /storage/software/Dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/storage/builab/Thibault/20241216_TetraCHE12over_TS/Doublet_STA/';


%%%%%%% Variables subject to change %%%%%%%%%%%
pixelSize = 14.00;
boxSize = 64;
filamentListFile = 'filamentList.csv';
alnDir = sprintf('%sintraAln', prjPath);
particleDir = sprintf('%sparticles', prjPath);
previewDir =[alnDir '/preview']; % created previously
mw = 12; % Number of parallel workers to run
gpu = [0]; % Alignment using gpu
initRefFile = 'templates/doublet_8nm_14.00Apx.em';
coneFlip = 1; % Search for polarity. 1 is yes. Should be 1 default.
avgLowpass = 28; % Angstrom
alnLowpass = 28; % Angstrom
shiftLimit = [10 10 4]; % Limit XYZ in pixel. Z should be half of periodicity
newRefFile = 'average_intraAln.em';
filamentPolarityListFile = sprintf('%sfilamentListPolarity.csv', prjPath);
combinedPreviewsDir = fullfile(previewDir, 'combined_previews');

%%%%%%% Do not change anything under here %%%%%

filamentList = readcell(filamentListFile, 'Delimiter', ',');
noFilament = length(filamentList);
alnLowpassPix = round(pixelSize/alnLowpass*boxSize);

template = dread(initRefFile);
newTemplate = zeros(boxSize, boxSize, boxSize);

filamentPolarityList = {};

% Need to go into alnDir to read the intraAln project
cd(alnDir)

% Use a containers.Map to accumulate aligned particles for each base filament
combinedMap = containers.Map('KeyType','char','ValueType','any');
if ~isfolder(combinedPreviewsDir)
    mkdir(combinedPreviewsDir);
end
% Initialize variables for combined previews
previousBaseName = ''; % To track the model being processed
partialSum = zeros(boxSize, boxSize, boxSize); % Sum of aligned particles for the current model
count = 0; % Number of filaments processed for the current model

%% Calculate the alignment of the filamentAverage to the initial reference
% transform the corresponding table for all particles
for idx = 1:noFilament
%for idx = 1:10
	aPath = ddb([filamentList{idx} ':a']); % Read the path of the alignment project average
	tPath = ddb([filamentList{idx} ':rt']);
	filamentAvg = dread(aPath);
  	disp(['Align ' filamentList{idx}]);
  	sal = dalign(dynamo_bandpass(filamentAvg,[1 alnLowpassPix]), dynamo_bandpass(template,[1 alnLowpassPix]),'cr',10,'cs',5,'ir',360,'is',10,'dim',boxSize, 'limm',1,'lim',shiftLimit,'rf',2,'rff',2, 'cone_flip', coneFlip); % cone_flip
	
	% Write out the transform
	writematrix([sal.p_shifts sal.p_eulers], [particleDir '/' filamentList{idx} '/xform.tbl'], 'Delimiter', 'tab', 'FileType', 'text');
	
	% Write out preview
	newTemplate = newTemplate + sal.aligned_particle;
	filt_aligned_particle = dynamo_bandpass(sal.aligned_particle, [1 round(pixelSize/avgLowpass*boxSize)]);
	img = sum(filt_aligned_particle(:,:,floor(boxSize/2) - 10: floor(boxSize/2) + 10), 3);
	% .png preview is rotated 90 deg clockwise, so we rotate 90 deg CC
    %img_rotated = imrotate(img, 90);
    img_rotated = flipud(img');
    imwrite(mat2gray(img_rotated), [previewDir '/' filamentList{idx} '_aln.png']);
	% Read last table from alignment
	tFilament = dread(tPath);
	% Read last transformation & applied to table
	tFilament_ali = dynamo_table_rigid(tFilament, sal.Tp);
	
	% Write table
	dwrite(tFilament_ali, [particleDir '/' filamentList{idx} '/aligned.tbl']);
	% Write aligned intraAvg
	dwrite(sal.aligned_particle, [alnDir '/avg/' filamentList{idx} '_aln.em']);
	
	% Read polarity here
	if abs(sal.p_eulers(2)) > 90
		% Flipping polarity
		filamentPolarityList{idx, 2} = 1;
	else
		% Not flipping polarity
		filamentPolarityList{idx, 2} = 0;
	end
	
	% Write Polarity list
	filamentPolarityList{idx, 1} = filamentList{idx};
    disp([filamentList{idx} ' Polarity ' num2str(filamentPolarityList{idx, 2})]);

    % To generate basename for each model (cilia)
    tokens = regexp(filamentList{idx}, '^([A-Za-z0-9]+_[0-9]{3})_[0-9]{1,2}$', 'tokens', 'once'); % Retreives basename of model
    if ~isempty(tokens)
        baseName = tokens{1};  % e.g. 'CHE12over_001'
    else
        warning(['Filament name "' filamentList{idx} '" does not match the expected pattern of CHE12over_000_0. Using entire name as baseName.']);
        baseName = filamentList{idx};
    end

    % If all filaments from one model have been added, generates a preview
    if (~strcmp(baseName, previousBaseName) && ~isempty(previousBaseName))
        % Compute and save the combined preview for the previous model
        avgMap = partialSum / count;
        filt_avg_map = dynamo_bandpass(avgMap, [1 round(pixelSize / alnLowpass * boxSize)]);
        img_combined = sum(filt_avg_map(:, :, floor(boxSize / 2) - 10 : floor(boxSize / 2) + 10), 3);

        % Rotate the combined image 90 degrees counterclockwise manually
        img_combined_rotated = flipud(img_combined'); % Equivalent to a 90-degree counterclockwise rotation

        % Save the combined preview image
        imwrite(mat2gray(img_combined_rotated), [combinedPreviewsDir '/' previousBaseName '_aln.png']);
        disp(['Saved combined preview for model: ' previousBaseName]);

        % Reset accumulators for the new model
        partialSum = zeros(boxSize, boxSize, boxSize);
        count = 0;
    end
    % Accumulate aligned particles for each base filament
    partialSum = partialSum + sal.aligned_particle;
    count = count + 1;

    % Update previousBaseName to the current filament's baseName
    previousBaseName = baseName;
end

% Handles the last model's preview
if count > 0 && ~isempty(previousBaseName)
    % Compute and save the combined preview for the last model
    avgMap = partialSum / count;
    filt_avg_map = dynamo_bandpass(avgMap, [1 round(pixelSize / alnLowpass * boxSize)]);
    img_combined = sum(filt_avg_map(:, :, floor(boxSize / 2) - 10 : floor(boxSize / 2) + 10), 3);

    % Rotate the combined image 90 degrees counterclockwise manually
    img_combined_rotated = flipud(img_combined'); % Equivalent to a 90-degree counterclockwise rotation

    % Save the combined preview image
    imwrite(mat2gray(img_combined_rotated), [combinedPreviewsDir '/' previousBaseName '_aln.png']);
    disp(['Saved combined preview for model: ' previousBaseName]);
end

cd ..

%% Write new list file and calculate new average
writecell(filamentPolarityList, filamentPolarityListFile);

newTemplate = newTemplate/noFilament;
dwrite(newTemplate, newRefFile)

