%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align repick particles
% and transform all the alignment to an updated table.
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /storage/software/Dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/storage/builab/Thibault/20241216_TetraCHE12over_TS/Doublet_STA/';

%% Input
boxSize = 64;
pixelSize = 14.00;
filamentRepickListFile = 'filamentRepickList.csv';
particleDir = sprintf('%sparticles_repick', prjPath);
mw = 6; % Number of parallel workers to run
gpu = [0]; % Alignment using gpu
template_name = 'templates/doublet_8nm_14.00Apx.em';
tableFileName = 'merged_particles_doublet_8nm_repick.tbl'; % merged particles table all
starFileName = 'merged_particles_doublet_8nm_repick.star'; % star file name for merged particles
tableOutFileName = 'merged_particles_doublet_8nm_repick.tbl'; % merged particles table all
pAlnAll = 'pAlnRepickParticles_8nm';
refMask = 'mask_doublet.em';
finalLowpass = 28; % Now implemented using in Angstrom
alnLowpass = 28; % Now implemented using Angstrom
zshift_limit = 3; % Should be half the periodicity, 4-nm for tip CP, 8-nm for doublet
newRefFile = 'average_allRepickParticles_doublet_8nm.em';

%%
filamentList = readcell(filamentRepickListFile, 'Delimiter', ',');
noFilament = length(filamentList);


template = dread(template_name);
 
% Combine all the particles into one table
% create table array
targetFolder = {};
tableName ={};

for idx = 1:noFilament
	targetFolder{idx} = [particleDir '/' filamentList{idx}];
	tableName{idx} = [particleDir '/' filamentList{idx} '/aligned.tbl'];
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
dvput(pAlnAll,'ite', [2 1]);
dvput(pAlnAll,'dim', [boxSize boxSize]);
dvput(pAlnAll,'low', [round(pixelSize/alnLowpass*boxSize) round(pixelSize/alnLowpass*boxSize)]);
dvput(pAlnAll,'cr', [15 6]);
dvput(pAlnAll,'cs', [5 2]);
dvput(pAlnAll,'ir', [15 6]);
dvput(pAlnAll,'is', [5 2]);
dvput(pAlnAll,'rf', [5 5]);
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
tPath = ddb([pAlnAll ':t:ite=last']); % This makes convertion to Relion better
dwrite(dread(tPath), tableOutFileName);
dwrite(dynamo_bandpass(a,[1 round(pixelSize/finalLowpass*boxSize)]),newRefFile);
