%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to align repick particles
% and transform all the alignment to an updated table.
% dynamoDMT v0.2b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set column 23 helicalID = 1

%%%%%%%% Before Running Script %%%%%%%%%%
%% Activate Dynamo
%run /london/data0/software/dynamo/dynamo_activate.m
run /data2/apps/dynamo/1.1.546/dynamo_activate.m 

% Change path to the correct directory
prjPath = '/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/tipCP_STA/';


%% Input
boxSize = 60;
pixelSize = 14.00;
refID = 2; % 1 for C2, 2 for C2
filamentRepickListFile = 'filamentRepickCPListC2.csv';
particleDir = sprintf('%sparticles_repick', prjPath);
mw = 8; % Number of parallel workers to run
gpu = [0 1]; % Alignment using gpu, for lots of particles, use more than 1
template_name = 'templates/tipC2_14.00Apx_new.em';
tableFileName = 'merged_particles_tipC2_repick.tbl'; % merged particles table all
starFileName = 'merged_particles_tipC2_repick.star'; % star file name for merged particles
tableOutFileName = 'merged_particles_tipC2_repick_align.tbl'; % merged particles table all
pAlnAll = 'pAlnRepickParticles_C2';
refMask = 'mask_cylinder.em'; % Better be tipC2/tip specific mask
finalLowpass = 28; % Now implemented using in Angstrom
alnLowpass = 28; % Now implemented using Angstrom
zshift_limit = 3; % Should be half the periodicity, 4-nm for tip CP, 8-nm for doublet
newRefFile = 'average_allRepickParticles_tipC2.em';

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

% Set the helicalTubeID (column 23) = refID
tMergedClean(:, 23) = refID;

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
