% getDataIntr.m, Youngcho Kim, Narayanan's Lab 2014
% Interactively get MedPC data parsed with mpcParser.m
% Parses MedPC files of requested file, retrun a structure with all MPC data
% Usage: getDataIntr(foldername); 
%
% the parseMPC is the return structure with metadata(date, subject, MSN)
% and timing data (S,H, etc)

% History
% 1/9/14  : quickly parse for MSN only, ask what to analyze.
% 1/19/14 : allow to select multiple MSN
% 1/23/14 : add error messages

% Get data folder, then MPC file list
% find the matching MSN, put them in a structure

function mpcParsed = getDataIntr(folder_name, MSN)
% get MEDPC file list from current directory by matching !2013-01-13
tmp_fileList = dir(folder_name);
fileList = regexpi({tmp_fileList.name},'(^\d{4}-\d{2}-\d{2}.txt|^!\d{4}-\d{2}-\d{2}.txt)','match');
% fileList = regexpi({tmp_fileList.name},'^d{4}-\d{2}-\d{2}','match');
fileList = [fileList{:}];

if size(fileList,2)==0
    error('No Matching MEDPC DATA file found!')
end

% only parse requested files, if group is defined
% if Date is 'ALL', need to parse every file
% if Date is defined, check file exist, then parse only requested files
% if group is available, MSN should be defined as well.


% quickly parse MSN only from file list
msnListNum =1;
msnListByFile= {};
for filenum=1:length(fileList)
    fString = fileread(fullfile(folder_name,char(fileList(filenum)))); % read the file as one string
    msns = regexp(fString,'MSN:\s([\w-]+)\r','tokens');
    msnListByFileNum = 1;
    for msnNum = 1:length(msns)
        msnListByFile(filenum,msnListByFileNum)=msns{msnNum};
        msnList(msnListNum) = msns{msnNum};
        msnListNum = msnListNum + 1;
        msnListByFileNum = msnListByFileNum + 1;
    end
end
if size(msnListByFile,1)==0
    error('No MEDPC DATA found from files!')
end

if ~exist('MSN','var') || isempty(MSN)
    % get all experiment protocol (MSN) from data, ask what to analyze
    % sort by frequency before ask 
    [uniqueMSN, ~, j] = unique(msnList);
    % if only one MSN is available, don't ask
    if size(uniqueMSN,2) == 1
        MSN = uniqueMSN(1);
    else
        [~, sorted_locations] = sort(hist(j, 1:max(j)));
        uniqueMSN = uniqueMSN(fliplr(sorted_locations)).';
        [Selection,~] = listdlg('PromptString','Select a Protocol', 'ListString',uniqueMSN);
        MSN = uniqueMSN(Selection);
    end
end
% find file list that have requested MSN
for msnSelected = 1:length(MSN)
    for filenum=1:length(fileList)
        tmp_msnbyFile = msnListByFile(filenum,:);
        tmp_msnbyFile = tmp_msnbyFile(~cellfun('isempty',tmp_msnbyFile));
        parsingFiles(msnSelected,filenum) = ismember(MSN{msnSelected},tmp_msnbyFile);
    end
end

% parse mpc files for requested MSN
mpcParsed = struct;
pfileList = fileList(any(parsingFiles,1));
for filenum=1:length(pfileList)
    fNAME = fullfile(folder_name,char(pfileList(filenum)));
    mpcParsed = mpcParser(fNAME,MSN,mpcParsed);
end

% % remove fields with no data 
% parsedFields = fieldnames(mpcParsed);
% for fieldNum = 1:length(parsedFields)
%     emptyData = 0;
%     for dataNum = 1:size(mpcParsed,2)
%         emptyData = emptyData + isempty(mpcParsed(dataNum).(parsedFields{fieldNum}));
%     end
%     if emptyData == size(mpcParsed,2)
%         mpcParsed = rmfield(mpcParsed,parsedFields{fieldNum});
%     end
% end
% 
% % data check for repeated testing in a day
% uniqueSbj = unique({mpcParsed.Subject});
% for i = 1:length(uniqueSbj)
%    subIdx = strcmp(uniqueSbj(i),{mpcParsed.Subject});
%    dateList = {mpcParsed(subIdx).StartDate};
%    if size(dateList,2) > size(unique(dateList),2)
%        warning('duplicated Animal ID %s found, check data!', char(uniqueSbj(i)));
%    end
% end

% remove data with no subject name
mpcParsed = mpcParsed(~strcmp('0',{mpcParsed.Subject}));

% fix Animal ID
for i=1:size(mpcParsed,2)
    sbj = mpcParsed(i).Subject;
    % check start with letter
    if ~isstrprop(sbj(1), 'alpha')
        sbj = strcat('FIX',sbj);
    end
    % replace "-" with "_"
    sbj = regexprep(sbj,'-','_');
    mpcParsed(i).Subject = sbj;
end

% Sort structure by MSN, subject name and timestamp(start date)
[~, ind] = sortrows([{mpcParsed.MSN}',{mpcParsed.Subject}',{mpcParsed.timeStamp}'], [1 2 3]);
mpcParsed = mpcParsed(ind);

