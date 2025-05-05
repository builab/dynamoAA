%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align all particles with a reference
% dynamoDMT v0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Before Running Script %%%%%%%%%%
%% Activate Dynamo
run /data2/apps/dynamo/1.1.546/dynamo_activate.m
%run /storage/software/Dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/mnt/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/Singlet_STA/';


%%%%%%% Variables subject to change %%%%%%%%%%%
pixelSize = 14.00;
boxSize = 60;
filamentListFile = 'filamentListPolarity.csv'; % Use manualFilamentCPList.csv if manual fixed, otherwise, filamentCPList.csv
alnDir = sprintf('%sintraAln', prjPath);
particleDir = sprintf('%sparticles', prjPath);
mw = 12; % Number of parallel workers to run
gpu = [0 1 2 3]; % Alignment using gpu
template_name = 'average_intraAvg_singlet.em'; % If you have a better reference, use it instead
tableFileName = 'merged_particles_singlet.tbl'; % merged particles table all
tableOutFileName = 'merged_particles_singlet_align.tbl'; % merged particles table all
starFileName = 'merged_particles_singlet.star'; % star file name for merged particles
pAlnAll = 'pAlnAllParticles';
refMask = 'mask_cylinder.em';
finalLowpass = 28; % Now implemented using in Angstrom
alnLowpass = 30; % Now implemented using Angstrom
zshift_limit = 4; % 4nm shift limit
outputFile = 'average_allParticles_singlet.em';


%%%%%%% Do not change anything under here %%%%%
filamentList = readcell(filamentListFile, 'Delimiter', ',');
noFilament = length(filamentList);

template = dread(template_name);
 
% Combine all the particles into one table
% create table array
targetFolder = {};
tableName ={};

for idx = 1:noFilament
	targetFolder{idx} = [particleDir '/' filamentList{idx, 1}];
	tableName{idx} = [particleDir '/' filamentList{idx, 1} '/aligned.tbl'];
end


% create ParticleListFile object (this object only exists temporarily in matlab)
plfClean = dpkdata.containers.ParticleListFile.mergeDataFolders(targetFolder,'tables',tableName);

% create and write the .star file
plfClean.writeFile(starFileName)

% create merged table
tMergedClean = plfClean.metadata.table.getClassicalTable();
dwrite(tMergedClean,tableFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform alignment of all particles with the ref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dcp.new(pAlnAll,'t', tableFileName, 'd', targetFolder{1}, 'template', template_name, 'masks','default','show',0, 'forceOverwrite',1);
dvput(pAlnAll,'data',starFileName)
dvput(pAlnAll,'file_mask',refMask)

% set alignment parameters
dvput(pAlnAll,'ite', [1 1]);
dvput(pAlnAll,'dim', [boxSize boxSize]); % Integer division of box size
dvput(pAlnAll,'low', [round(boxSize*pixelSize/alnLowpass) round(boxSize*pixelSize/alnLowpass)]); % Low pass filter
dvput(pAlnAll,'cr', [15 6]);
dvput(pAlnAll,'cs', [5 2]);
dvput(pAlnAll,'ir', [15 6]);
dvput(pAlnAll,'is', [5 2]);
dvput(pAlnAll,'rf', [2 2]);
dvput(pAlnAll,'rff', [2 2]);
dvput(pAlnAll,'lim', [zshift_limit zshift_limit]);
dvput(pAlnAll,'limm',[1 2]);
dvput(pAlnAll,'sym', 'c1'); % 
    
% set computational parameters
dvput(pAlnAll,'dst','matlab_gpu','cores',1,'mwa',mw);
dvput(pAlnAll,'gpus',gpu);

% check/unfold/run
dvrun(pAlnAll,'check',true,'unfold',true);

aPath = ddb([pAlnAll ':a']);
a = dread(aPath);
tPath = ddb([pAlnAll ':t:ite=last']); % This is correct but might not be prone to more error!!!
dwrite(dread(tPath), tableOutFileName);
dwrite(dynamo_bandpass(a,[1 round(boxSize*pixelSize/finalLowpass)]), outputFile);
