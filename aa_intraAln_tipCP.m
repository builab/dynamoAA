%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align subtomogram within the same filament using twist
% dynamoDMT v1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Before Running Script %%%%%%%%%%%%%%%
%%% Activate Dynamo
run /data2/apps/dynamo/1.1.546/dynamo_activate.m
%run /storage/software/Dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/storage/builab/Thibault/20241216_TetraCHE12over_TS/TipCP_STA/';

%%%%%%% Variables subject to change %%%%%%%%%%%

docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath);
filamentListFile = sprintf('%sfilamentList.csv', prjPath);
alnDir = sprintf('%sintraAln', prjPath);
particleDir = sprintf('%sparticles', prjPath);
boxSize = 60; % Original extracted subvolume size
mw = 12; % Number of parallel workers to run
gpu = [0]; % Alignment using gpu for titann setting
pixelSize = 14.00; % Angstrom per pixel
avgLowpass = 30; % In Angstrom to convert to Fourier Pixel
alnLowpass = 30; % In Angstrom to convert to Fourier Pixel, 50 if you only need to align the filament well, 35-40Angstrom for clear tubulin MT
zshift_limit = 4; % ~4nm 

%%%%%%% Do not change anything under here %%%%%

%% Generate an initial reference average for each filament
filamentList = readcell(filamentListFile, 'Delimiter', ',');

mkdir(alnDir);
mkdir([alnDir '/avg']); %filter averages
mkdir([alnDir '/preview']); % preview images

%% Loop through filamentList
cd(alnDir)
for idx = 1:length(filamentList)
    tableName = [particleDir '/' filamentList{idx} '/crop.tbl'];
    tOri = dread(tableName);
    template = [particleDir '/' filamentList{idx} '/template.em'];
    prjPaticlesDir = [particleDir '/' filamentList{idx}];
    prj_intra = [filamentList{idx}];    

    % create alignment project
    dcp.new(prj_intra,'d',prjPaticlesDir,'t',tableName, 'template', template, 'masks','default','show',0);

    % set alignment parameters for 2 rounds
    dvput(prj_intra,'ite', [2]); % no iterations 3 is reasonable
    dvput(prj_intra,'dim', [boxSize]); % Use 1/2 box size for quicker but full size for good res
    dvput(prj_intra,'low', [round(pixelSize/alnLowpass*boxSize)]); % lowpass filter
    dvput(prj_intra,'cr', [15]); % cone range
    dvput(prj_intra,'cs', [5]); % cone search step
    dvput(prj_intra,'ir', [15]); % inplane rotation
    dvput(prj_intra,'is', [5]); % inplane search step
    dvput(prj_intra,'rf', [5]); % refinement
    dvput(prj_intra,'rff', [2]); % refinement factor
    dvput(prj_intra,'lim', [zshift_limit]); % shift limit
    dvput(prj_intra,'limm',[1]); % limit mode
    dvput(prj_intra,'sym', 'c1'); % symmetry

    % set computational parameters
    dvput(prj_intra,'dst','matlab_gpu','cores',1,'mwa',mw);
    dvput(prj_intra,'gpus',gpu);
    
    %CPU
    %dvput(prj_intra,'dst', 'matlab_parfor','cores',12,'mwa',mw);

    % check/unfold/run
    dvrun(prj_intra,'check',true,'unfold',true);
    
    % Generate the average & filter to 30 Angstrom & a png preview
    aPath = ddb([filamentList{idx} ':a']); % Read the path of the alignment project average
    filamentAvg = dread(aPath);
    filamentAvg = dynamo_bandpass(filamentAvg,[1 round(pixelSize/avgLowpass*boxSize)]);
    dwrite(filamentAvg, ['avg/' filamentList{idx} '.em']);
    
    % Preview
    img = sum(filamentAvg(:,:,floor(boxSize/2) - 5: floor(boxSize/2) + 5), 3);
    % .png preview is rotated 90 deg clockwise, so we rotate 90 deg CC
    img_rotated = flipud(img');
    imwrite(mat2gray(img_rotated), ['preview/' filamentList{idx} '.png'])
end

cd ..
