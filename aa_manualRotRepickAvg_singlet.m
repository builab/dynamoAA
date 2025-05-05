%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to manual edit the xform file to rotate the filament
% Auto edit the rotation angle in the 6 column of the xform file
% -angle from imod mean + angle here as no polarity is inverted
% DynamoMT v0.1 (identical for 13 & 14 PF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%% Before Running Script %%%%%%%%%%
%%% Activate Dynamo
run /data2/apps/dynamo/1.1.546/dynamo_activate.m
%run /storage/software/Dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/mnt/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/Singlet_STA/';


%%%%%%% Variables subject to change %%%%%%%%%%%

pixelSize = 14.00;
boxSize = 60;
filamentListFileManualRot= 'filamentRepickListManualRot.csv';
alnDir = sprintf('%sintraAln_repick', prjPath);
particleDir = sprintf('%sparticles_repick', prjPath);
previewDir =[alnDir '/preview']; % created from previously
mw = 2; % Number of parallel workers to run
avgLowpass = 28; % Angstrom
skipIntraAln = 0; % In repick, you can opt to use intraRepickAln or not
useMidRegionOnly = 0; % Use only the middle section to preserve features


%%%%%%% Do not change anything under here %%%%%

filamentList = readcell(filamentListFileManualRot, 'Delimiter', ',');
noFilament = size(filamentList, 1);

% Need to go into alnDir to read the intraAln project
cd(alnDir)

%% Calculate the alignment of the filamentAverage to the initial reference
% transform the corresponding table for all particles
for idx = 1:noFilament
	if skipIntraAln > 0
		aPath = ([particleDir '/' filamentList{idx, 1} '/template.em']); % Read the path of the alignment project average
		tPath = ([particleDir '/' filamentList{idx, 1} '/crop.tbl']); 
    	else
        	if useMidRegionOnly > 0
            		filamentAvg = dread([alnDir '/avg/' filamentList{idx, 1} '_mid.em']);
        	else
            		filamentAvg = dread(ddb([filamentList{idx, 1} ':a'])); % Read the path of the alignment project average
        	end
        	tPath = ddb([filamentList{idx, 1} ':rt']);
	end
    	disp(filamentList{idx})

	t_xform = load([particleDir '/' filamentList{idx,1} '/xform.tbl']);
	t_xform(1, 6) = t_xform(1, 6) - filamentList{idx,2}; % Reverse the angle from filamentListManualRot.csv
	Tp.type = 'shiftrot';
	Tp.shifts = t_xform(1, 1:3);
	Tp.eulers = t_xform(1, 4:6);
	
	% Read last table from alignment
	tFilament = dread(tPath);
	% Read last transformation & applied to table
	tFilament_ali = dynamo_table_rigid(tFilament, Tp);
	% OverWrite the old table
	dwrite(tFilament_ali, [particleDir '/' filamentList{idx,1} '/aligned.tbl']);
	% Write aligned intraAvg
	average = dpkgeom.rotation.smoothShiftRot(filamentAvg, Tp.shifts, Tp.eulers);
	dwrite(dynamo_bandpass(average, [1 round(pixelSize/avgLowpass*boxSize)]), [alnDir '/avg/' filamentList{idx,1} '_manual_aln.em']);

end
 
cd ..


