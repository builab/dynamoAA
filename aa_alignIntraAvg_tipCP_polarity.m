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
prjPath = '/mnt/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/tipCP_STA/';


%%%%%%% Variables subject to change %%%%%%%%%%%
pixelSize = 14.00;
boxSize = 60;
filamentListFile = 'filamentList.csv';
alnDir = sprintf('%sintraAln', prjPath);
particleDir = sprintf('%sparticles', prjPath);
previewDir =[alnDir '/preview']; % created previously
mw = 12; % Number of parallel workers to run
gpu = [0]; % Alignment using gpu
initRefFile = 'templates/tipC1_14.00Apx_new.em';
coneFlip = 0; % Search for polarity. 1 is yes. Should be 0 default for this script
avgLowpass = 28; % Angstrom
alnLowpass = 28; % Angstrom
shiftLimit = [10 10 4]; % Limit XYZ in pixel. Z should be half of periodicity
newRefFile = 'average_intraAln.em';
filamentPolarityListFile = sprintf('%sfilamentListPolarity.csv', prjPath);
maskFile = sprintf('%smask_cylinder.em', prjPath); % Cylinder mask, must be quite big
ciliaPolarityFile = sprintf('%sciliaPolarity.csv', prjPath);

%%%%%%% Do not change anything under here %%%%%

filamentList = readcell(filamentListFile, 'Delimiter', ',');
noFilament = length(filamentList);
ciliaPolarity = readcell(ciliaPolarityFile);
alnLowpassPix = round(pixelSize/alnLowpass*boxSize);

template = dread(initRefFile);
newTemplate = zeros(boxSize, boxSize, boxSize);

filamentPolarityList = {};

% Need to go into alnDir to read the intraAln project
cd(alnDir)

%% Calculate the alignment of the filamentAverage to the initial reference
% transform the corresponding table for all particles
for idx = 1:noFilament
	aPath = ddb([filamentList{idx} ':a']); % Read the path of the alignment project average
	tPath = ddb([filamentList{idx} ':rt']);
	filamentAvg = dread(aPath);
	
	
	% Get the tomoName (baseName)
    tokens = regexp(filamentList{idx}, '^(.*)_[0-9]{1,2}$', 'tokens', 'once'); % Retreives basename of model
    if ~isempty(tokens)
        baseName = tokens{1};  % e.g. 'CHE12over_001'
    else
        warning(['Filament name "' filamentList{idx} '" does not match the expected pattern of XXX_YYY_D. Using entire name as baseName.']);
        baseName = filamentList{idx};
    end
	
	% Check the polarity 
	matchIndex = strcmp(ciliaPolarity(:,1), baseName);
	
	% Get the corresponding polarity value from the second column
	if any(matchIndex)
    	polarity = ciliaPolarity{matchIndex, 2};
    	disp(['Polarity for ', baseName, ' is: ', num2str(polarity)]);
	else
    	disp('String not found in ciliaPolarity.');
	end

	% If polarity = 1, multiply with [0 180 0], then align
	Tflip.type = 'shiftrot';
	Tflip.shifts = [0 0 0];
	Tflip.eulers = [0 180 0];
	if polarity > 0
		Tflip.euler = [0 180 0];
		flippedVol = dpkgeom.rotation.smoothShiftRot(filamentAvg, Tflip.shifts, Tflip.eulers);
	else
		flippedVol = filamentAvg;
	end
	
  	disp(['Align ' filamentList{idx}]);
  	sal = dalign(dynamo_bandpass(flippedVol,[1 alnLowpassPix]), dynamo_bandpass(template,[1 alnLowpassPix]),'mask', maskFile, 'cr',10,'cs',5,'ir',360,'is',10,'dim',boxSize, 'limm',1,'lim',shiftLimit,'rf',2,'rff',2, 'cone_flip', coneFlip); % cone_flip
	
	% Write out the transform
	writematrix([sal.p_shifts sal.p_eulers], [particleDir '/' filamentList{idx} '/xform.tbl'], 'Delimiter', 'tab', 'FileType', 'text');
	
	
	% Write out preview
	newTemplate = newTemplate + sal.aligned_particle;
	filt_aligned_particle = dynamo_bandpass(sal.aligned_particle, [1 round(pixelSize/avgLowpass*boxSize)]);
	img = sum(filt_aligned_particle(:,:,floor(boxSize/2) - 10: floor(boxSize/2) + 10), 3);
	% .png preview is rotated 90 deg clockwise, so we rotate 90 deg CC
    img_rotated = flipud(img');
    imwrite(mat2gray(img_rotated), [previewDir '/' filamentList{idx} '_aln.png']);
    
	% Read last table from alignment
	tFilament = dread(tPath);
	% Read last transformation & applied to table, take care of the flip as well
	tFilament_flip = dynamo_table_rigid(tFilament, Tflip)
	tFilament_ali = dynamo_table_rigid(tFilament_flip, sal.Tp);
	
	% Write table
	dwrite(tFilament_ali, [particleDir '/' filamentList{idx} '/aligned.tbl']);
	% Write aligned intraAvg
	dwrite(sal.aligned_particle, [alnDir '/avg/' filamentList{idx} '_aln.em']);
	
	% Write polarity here (for compatibility)
	filamentPolarityList{idx, 2} = polarity;
	
	% Write Polarity list
	filamentPolarityList{idx, 1} = filamentList{idx};
    disp([filamentList{idx} ' Polarity ' num2str(filamentPolarityList{idx, 2})]);


end


cd ..

%% Write new list file and calculate new average
writecell(filamentPolarityList, filamentPolarityListFile);

newTemplate = newTemplate/noFilament;
dwrite(newTemplate, newRefFile)


