%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align filament to a reference
% dynamoMT v0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The new list will contain 2 columns (Filament, polarity)

%%%%%%%% Before Running Script %%%%%%%%%%%%%%%

%%% Activate Dynamo
run /data2/apps/dynamo/1.1.546/dynamo_activate.m
%run /storage/software/Dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/mnt/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/Singlet_STA/';


%%%%%%% Variables subject to change %%%%%%%%%%%
pixelSize = 14.00;
boxSize = 60;
filamentListFile = 'filamentRepickList.csv';
alnDir = sprintf('%sintraAln_repick', prjPath);
particleDir = sprintf('%sparticles_repick', prjPath);
previewDir =[alnDir '/preview']; % created previously
mw = 12; % Number of parallel workers to run
gpu = [0 1]; % Alignment using gpu
initRefFile = 'templates/Singlet_SPEF1_14.00Apx.em';
coneFlip = 0; % Search for polarity. 1 is yes. Should be 0 default.
avgLowpass = 28; % Angstrom
alnLowpass = 28; % Angstrom
shiftLimit = [3 3 3]; % Limit XYZ in pixel. Z should be half of periodicity
newRefFile = 'average_RepickAvg.em';
maskFile = sprintf('%mask_singlet.em', prjPath);
combinedPreviewsDir = fullfile(previewDir, 'combined_previews');

%%%%%%% Do not change anything under here %%%%%

filamentList = readcell(filamentListFile, 'Delimiter', ',');
noFilament = length(filamentList);
alnLowpassPix = round(pixelSize/alnLowpass*boxSize);

template = dread(initRefFile);
newTemplate = zeros(boxSize, boxSize, boxSize);

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
	aPath = ddb([filamentList{idx} ':a']); % Read the path of the alignment project average
	tPath = ddb([filamentList{idx} ':rt']);
	filamentAvg = dread(aPath);
  	disp(['Align ' filamentList{idx}]);
  	sal = dalign(dynamo_bandpass(filamentAvg,[1 alnLowpassPix]), dynamo_bandpass(template,[1 alnLowpassPix]),'mask', maskFile, 'cr',10,'cs',5,'ir',360,'is',10,'dim',boxSize, 'limm',1,'lim',shiftLimit,'rf',2,'rff',2, 'cone_flip', coneFlip); % cone_flip
	
	% Write out the transform
	writematrix([sal.p_shifts sal.p_eulers], [particleDir '/' filamentList{idx} '/xform.tbl'], 'Delimiter', 'tab', 'FileType', 'text');
	
	% Write out preview
	newTemplate = newTemplate + sal.aligned_particle;
	filt_aligned_particle = dynamo_bandpass(sal.aligned_particle, [1 round(pixelSize/avgLowpass*boxSize)]);
	img = sum(filt_aligned_particle(:,:,floor(boxSize/2) - 10: floor(boxSize/2) + 10), 3);
    	img_rotated = flipud(img'); % Rotate CC 90 deg
	imwrite(mat2gray(img_rotated), [previewDir '/' filamentList{idx} '_aln.png']);
	% Read last table from alignment
	tFilament = dread(tPath);
	% Read last transformation & applied to table
	tFilament_ali = dynamo_table_rigid(tFilament, sal.Tp);
	
	% Write table
	dwrite(tFilament_ali, [particleDir '/' filamentList{idx} '/aligned.tbl']);
	% Write aligned intraAvg
	dwrite(sal.aligned_particle, [alnDir '/avg/' filamentList{idx} '_aln.em']);
	
    % To generate basename for each model (cilia)
    tokens = regexp(filamentList{idx}, '^(.*)_[0-9]{1,2}$', 'tokens', 'once'); % Retreives basename of model
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
newTemplate = newTemplate/noFilament;
dwrite(newTemplate, newRefFile)

