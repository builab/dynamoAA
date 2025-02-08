%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the project
% dynamoDMT v1.0 (now working with Warp workflow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import tomogram & create catalog using the GUI
% required IMOD
% Need to copy all the .tlt file to the same folder as the reconstructions


%%%%%% Instructions %%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Change the following paths to the right directories.
%%%% modDir contains the files with the tomograms and path references the root directory for analysis
% Step 2: Indicate the delimiter that references all mod files of interest.
%%%% CU428*/tip_CP.mod
% Step 3: Indicate the suffix that needs to be removed in the various .mod files.
%%%% In the same way, indicate the suffix that needs to be added to reference the mrc files.
% Step 4: Change the pixel size
%%%% 

%%%%%% Activate Dynamo %%%%%%%%%%%%%%%%%%%%%%%%
run /storage/software/Dynamo/dynamo_activate.m

%%%%%%% Variables subject to change %%%%%%%%%%%

modDir = '/storage/builab/Thibault/20241216_TetraCHE12over_TS/warp_tiltseries/reconstruction/';
prjPath = '/storage/builab/Thibault/20241216_TetraCHE12over_TS/Doublet_STA/';
modFileDelimiter = '*_doublets.mod';
stringToBeRemoved = '_14.00Apx_doublets.mod';
mrcFileDelimiter = '"*_14.00Apx.mrc"';
mrcStringToBeRemoved = '_14.00Apx.mrc';
recSuffix = '_14.00Apx';
apixel = '14.00';
useWarpRec = 1;

%%%%%%% Do not change anything under here %%%%%

pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'createModTxt.sh');
cmdStr = [pathToModelScript ' ' modDir ' ' modFileDelimiter ' ' prjPath];
system(cmdStr);
pathToMrcScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'createMrcTxt.sh');
cmdStr = [pathToMrcScript ' ' modDir ' ' mrcFileDelimiter ' ' prjPath];
system(cmdStr);

modelfile = sprintf('%smodfiles.txt', prjPath);
mrcfile = sprintf('%smrcfiles.txt', prjPath);
catalogs = sprintf('%scatalogs', prjPath);
listOfTomograms = sprintf('%scatalogs/listOfTomograms', prjPath);
c001Dir = sprintf('%scatalogs/c001', prjPath);
modelDir = sprintf('%smodels/', prjPath);
mkdir(catalogs);
mkdir(listOfTomograms);
mkdir(c001Dir);
mkdir(modelDir);

catPath = sprintf('%scatalogs/c001', prjPath);
docFilePath = sprintf('%scatalogs/tomograms.doc', prjPath); % This has to be consistent with the vll after catalog creation
vllFilePath = sprintf('%scatalogs/tomograms.vll', prjPath);

%% Populate the necessary files into a single location. Use this directory
% to create .vll and .doc files

% $./vllAndDocScript.sh modDir modelDestination modelfile stringTobeRemoved docFilePath vllFilePath 

if useWarpRec > 0
    pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'vllAndDocScriptWarp.sh');
else
	pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'vllAndDocScript.sh');
end
cmdStr = [pathToModelScript ' ' modDir ' ' listOfTomograms ' ' modelfile ' ' stringToBeRemoved ' ' docFilePath ' ' vllFilePath ' ' recSuffix '.mrc ' apixel];
system(cmdStr);

% Create new catalogue from vll file, delete old one
dcm('c', catPath, 'fromvll', vllFilePath, 'delete_old', 1)

% Create text files using model2point
% Imod coordinate should be in text file, clicking along the filament (no direction needed)
% model2point -Contour imodModel.mod imodModel.txt

pathToModelScript = fullfile(sprintf('%sdynamoDMT', prjPath), 'model2pointscript.sh');

cmdStr = [pathToModelScript ' ' modDir ' ' modelDir ' ' modelfile ' ' stringToBeRemoved];
system(cmdStr);
