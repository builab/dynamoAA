%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to sort MT into C1 or C2 PF with a polarity check
% dynamoAA v0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polarity check & compare to C1 and C2 ref 
% derived from dynamoMT/mt_sortPFs.m
% The new list will contain 3 columns (Filament, polarity, C1 or C2)
% Also, the 

%%%%%%%% Before Running Script %%%%%%%%%%%%%%%

%%% Activate Dynamo
run /data2/apps/dynamo/1.1.546/dynamo_activate.m
%run /storage/software/Dynamo/dynamo_activate.m

% Change path to the correct directory
prjPath = '/mnt/lima/huy/data0/20221128_TetraCU428Membrane_26k_TS/tipCP_STA/';


%%%%%%% Variables subject to change %%%%%%%%%%%
pixelSize = 14.00;
boxSize = 60;
filamentListFile = 'filamentRepickList.csv';
alnDir = sprintf('%sintraAln_repick', prjPath);
particleDir = sprintf('%sparticles_repick', prjPath);
previewDir =[alnDir '/preview']; % created previously
mw = 12; % Number of parallel workers to run
gpu = [0]; % Alignment using gpu
initRefFiles = {'templates/tipC1_14.00Apx_new.em', 'templates/tipC2_14.00Apx_new.em'};
refPFs = [1 2];
coneFlip = 0; % Search for polarity. Should be 0 here.
avgLowpass = 25; % Angstrom
alnLowpass = 25; % Angstrom
shiftLimit = [5 5 4]; % Limit XYZ in pixel. Z should be half of periodicity
newRefFile = 'average_sortRepickCP.em';
filamentCPListFile = sprintf('%sfilamentRepickCPList.csv', prjPath);
CClistFile = sprintf('%sCCListRepick.csv', prjPath);
maskFile = sprintf('%smask_cylinder.em', prjPath); % Cylinder mask, must be quite big


%%%%%%% Do not change anything under here %%%%%

filamentList = readcell(filamentListFile, 'Delimiter', ',');
noFilament = length(filamentList);
alnLowpassPix = round(pixelSize/alnLowpass*boxSize);
newTemplate = {};

template = {};
for refIdx = 1:length(initRefFiles)
	template{refIdx} = dread(initRefFiles{refIdx});
	newTemplate{refIdx} = zeros(boxSize, boxSize, boxSize);
end

filamentCPList = {};


% Need to go into alnDir to read the intraAln project
cd(alnDir)

%% Calculate the alignment of the filamentAverage to the initial reference
% transform the corresponding table for all particles
% For CClist, we will have tomoNo,filamentNo, CC_ref1, CC_ref2

CClist = zeros(noFilament,  5 + length(template));

for idx = 1:noFilament
	aPath = ddb([filamentList{idx} ':a']); % Read the path of the alignment project average
	tPath = ddb([filamentList{idx} ':rt']);
	tFilament = dread(tPath);
	filamentAvg = dread(aPath);
  	disp(['Align ' filamentList{idx}]);
  	sal = {};
  	maxCC = 0;
  	maxPF = 0;
	CClist(idx, 1:2) = tFilament(1, [20 23]);
  	for refIdx = 1:length(template)
  		sal{refIdx} = dalign(dynamo_bandpass(filamentAvg,[1 alnLowpassPix]), dynamo_bandpass(template{refIdx},[1 alnLowpassPix]),'mask', maskFile, 'cr', 9,'cs',3,'ir',15,'is',3,'dim',boxSize, 'limm',1,'lim',shiftLimit,'rf',2,'rff',2, 'cone_flip', coneFlip); % cone_flip
		
		if sal{refIdx}.ccmax(end) > maxCC
			maxCC = sal{refIdx}.ccmax(end);
			maxPF = refIdx;
		end
		CClist(idx, 2 + refIdx) = sal{refIdx}.ccmax(end);
	end
	
	% Write out the transform
	writematrix([sal{maxPF}.p_shifts sal{maxPF}.p_eulers], [particleDir '/' filamentList{idx} '/xform.tbl'], 'Delimiter', 'tab', 'FileType', 'text');
	
	% Write out preview
	newTemplate{maxPF} = newTemplate{maxPF} + sal{maxPF}.aligned_particle;
	filt_aligned_particle = dynamo_bandpass(sal{maxPF}.aligned_particle, [1 round(pixelSize/avgLowpass*boxSize)]);
	img = sum(filt_aligned_particle(:,:,floor(boxSize/2) - 10: floor(boxSize/2) + 10), 3);
	% .png preview is rotated 90 deg clockwise, so we rotate 90 deg CC
    	img_rotated = flipud(img');
	imwrite(mat2gray(img), [previewDir '/' filamentList{idx} '_aln_C' num2str(refPFs(maxPF)) '.png']);
	% Read last table from alignment
	%tFilament = dread(tPath);
	% Read last transformation & applied to table
	tFilament_ali = dynamo_table_rigid(tFilament, sal{maxPF}.Tp);
	% Write table
	dwrite(tFilament_ali, [particleDir '/' filamentList{idx} '/aligned.tbl']);
	% Write aligned intraAvg
	dwrite(sal{maxPF}.aligned_particle, [alnDir '/avg/' filamentList{idx} '_aln_C' num2str(refPFs(maxPF)) '.em']);
	
	% Read polarity here
	if abs(sal{maxPF}.p_eulers(2)) > 90
		% Flipping polarity
		filamentCPList{idx, 2} = 1;
	else
		% Not flipping polarity
		filamentCPList{idx, 2} = 0;
	end
	
	disp([filamentList{idx} ' Polarity ' num2str(filamentCPList{idx, 2}) ' and class C' num2str(refPFs(maxPF))]);
	% Write PF list
	filamentCPList{idx, 1} = filamentList{idx};
	filamentCPList{idx, 3} = refPFs(maxPF); 
	   
end
 
cd ..

% ReCheck the assignment
% Sample CClist matrix (n x 7), where columns are: tomoNo, filamentNo, CC1, CC2, CP, CC1-CC2, FinalCP
% CClist = [tomoNo, filamentNo, CC1, CC2, CP, abs(CC1-CC2), FinalCP];
% Compute the 5th column: 1 if CC1 > CC2, otherwise 2
CClist(:,5) = (CClist(:,4) > CClist(:,3)) + 1;  

% Compute the 6th column: max of CC1 and CC2
CClist(:,6) = abs(CClist(:,3) - CClist(:,4));

% Step 1: Initialize column 7 with the same values as column 5
CClist(:,7) = CClist(:,5);

% Step 2: Find unique tomoNo
uniqueTomo = unique(CClist(:,1));

% Step 3: Iterate over each tomoNo
for i = 1:length(uniqueTomo)
    tomoID = uniqueTomo(i);
    
    % Get indices of rows corresponding to this tomoNo
    idx = find(CClist(:,1) == tomoID);
    
    if length(idx) == 2  % Each tomoNo should have exactly 2 filamentNo
        % Extract rows
        row1 = idx(1);
        row2 = idx(2);
        
        % Step 4: Check if both rows have the same value in column 5
        if CClist(row1,5) == CClist(row2,5)
            % Find the row with the larger column 6 value, larger classification difference
            if CClist(row1,6) > CClist(row2,6)
                CClist(row2,7) = 3 - CClist(row2,5); % Swap 1 ↔ 2
            else
                CClist(row1,7) = 3 - CClist(row1,5); % Swap 1 ↔ 2
            end
        end
    end
end

writematrix(CClist, CClistFile);

% Output the final filament assignment
for idx = 1:noFilament
	filamentCPList{idx, 3} = CClist(idx, 7);
	%filamentCPList{idx, 4} = CClist(idx, 5);
end

%% Calculate average
writecell(filamentCPList, filamentCPListFile);

% Write separate list files for different PFs
numPF = cell2mat(filamentCPList(:, 3));
for refIdx = 1:length(refPFs)
	subFilamentList = filamentCPList(numPF == refPFs(refIdx), :);
	if isempty(subFilamentList)
		disp(['No filament with C' num2str(refPFs(refIdx))]);
	else
		writecell(subFilamentList, strrep(filamentCPListFile, '.csv', ['C' num2str(refPFs(refIdx)) '.csv']));
		newTemplate{refIdx} = newTemplate{refIdx}/length(subFilamentList);
		dwrite(newTemplate{refIdx}, strrep(newRefFile, '.em', ['_C' num2str(refPFs(refIdx)) '.em']))
	end
end
