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
prjPath = '/storage/builab/Thibault/20241216_TetraCHE12over_TS/tipCP_STA/';

%%%%%%% Variables subject to change %%%%%%%%%%%

pixelSize = 14.00;
boxSize = 60;
filamentListFileManualRot= 'filamentListPolarity.csv';
alnDir = sprintf('%sintraAln', prjPath);
particleDir = sprintf('%sparticles', prjPath);
previewDir =[alnDir '/preview']; % created from previously
mw = 2; % Number of parallel workers to run
avgLowpass = 28; % Angstrom
manualFilamentPolarityListFile = sprintf('%smanualFilamentCPList.csv', prjPath);


%%%%%%% Do not change anything under here %%%%%

filamentList = readcell(filamentListFileManualRot, 'Delimiter', ',');
noFilament = size(filamentList, 1);

filamentPolarityList = {};

% Need to go into alnDir to read the intraAln project
cd(alnDir)

%% Calculate the alignment of the filamentAverage to the initial reference
% transform the corresponding table for all particles
for idx = 1:noFilament
    filamentAvg = dread([alnDir '/avg/' filamentList{idx, 1} '_aln.em']); % Read the path of the alignment project average
    if ~exist([particleDir '/' filamentList{idx,1} '/manual_xform.tbl'], 'file') %Check is manual_xform.tbl file exists
        filamentPolarityList{idx, 1} = filamentList{idx, 1};
        filamentPolarityList{idx, 2} = filamentList{idx, 2};
        continue;
    end
	t_xform = load([particleDir '/' filamentList{idx,1} '/manual_xform.tbl']);
    Tp.type = 'shiftrot';
	Tp.shifts = t_xform(1, 1:3);
	Tp.eulers = t_xform(1, 4:6);

	% Read last table from alignment
    tPath = ([particleDir '/' filamentList{idx, 1} '/aligned.tbl']); 
	tFilament = dread(tPath);
	% Read manual transformation & applied to table and overwrite the
    % aligned.tbl
	tFilament_ali = dynamo_table_rigid(tFilament, Tp);
    dwrite(tFilament_ali, [particleDir '/' filamentList{idx,1} '/aligned.tbl']);

	% Read last table from alignment
	%vaverage = dpkgeom.rotation.smoothShiftRot(filamentAvg, Tp.shifts, Tp.eulers);
    vol = dpkgeom.rotation.smoothShiftRot(filamentAvg, Tp.shifts, Tp.eulers);
	% Write out
	dwrite(dynamo_bandpass(vol, [1 round(pixelSize/avgLowpass*boxSize)]), [alnDir '/avg/' filamentList{idx,1} '_manual_aln.em']);
    % Perhaps write out image for quick preview

    % Flips polarity if modified
	if t_xform(1,5) > 90
        if filamentList{idx,2} == 0
                filamentPolarityList{idx,2} = 1;
        else
            filamentPolarityList{idx,2} = 0;
        end
    else
        filamentPolarityList{idx,2} = filamentList{idx,2};
    end
    
    filamentPolarityList{idx, 1} = filamentList{idx, 1};

end
 
cd ..

writecell(filamentPolarityList, manualFilamentPolarityListFile);

