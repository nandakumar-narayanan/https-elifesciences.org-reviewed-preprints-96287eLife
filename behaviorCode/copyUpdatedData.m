function copyUpdatedData(data_source, local_data, MSN, sbjIDs, dateRange)
% data_source = 'D:\BehaviorDataSync\Data';
% local_data = 'C:\Users\youukim\Desktop\SimonTask\behaviorData';
% MSN={'SIMON_FR1'    'SIMON_FR1HOLD'    'SIMON_FR1HOLD_TO'};

if nargin < 5
    dateRange = [];
end
if nargin < 4
    sbjIDs = [];
end

if exist(local_data, 'dir') ~= 7
    fprintf('local data folder do not exit, create one\n')
    mkdir(local_data);
end

% find last date of data from local data folder
local_fileList = dir(local_data);
fileList = regexpi({local_fileList.name},'^!\d{4}-\d{2}-\d{2}','match','once');
fileList = fileList(~cellfun(@isempty, fileList));

% check whether currenct local data have matching 
if ~isempty(fileList)
    mfiles = getIncFiles(data_source, fileList, MSN, sbjIDs);
else
    mfiles = [];
end

% if no pre-existing data files or existing local data do not have matching data, check last 30 days
if ~isempty(fileList) && ~isempty(mfiles)
    dataDateStr = cellfun(@(x) x(2:end),fileList, 'UniformOutput', false);
    dataDate = datetime(dataDateStr,'InputFormat','yyyy-MM-dd');
    dataDate = sort(dataDate);
    localLastDate = dataDate(end);
else
    dataDate = datetime('today');
    localLastDate = dataDate - calmonths(1);
end

ind_dates = false;
if iscell(dateRange) && ~isempty(dateRange) % date range is requested
    if length(dateRange)<3 % 1 or 2 is considered as range
        localLastDate = datetime(dateRange{1},'InputFormat','yyyy-MM-dd');
        if length(dateRange)>1
            endDate = datetime(dateRange{2},'InputFormat','yyyy-MM-dd');
        elseif length(dateRange)==1
            endDate = datetime('today');
        else

        end
    else % more than 2 will be considered as individual dates
        ind_dates = true;
    end
else
    endDate = [];
end

% find updated data files from data source folder
source_fileList = dir(data_source);
fileList = regexpi({source_fileList.name},'^!\d{4}-\d{2}-\d{2}','match');
fileList = [fileList{:}];
dataDateStr = cellfun(@(x) x(2:end),fileList, 'UniformOutput', false);
dataDate = datetime(dataDateStr,'InputFormat','yyyy-MM-dd');
[dataDate, sortIdx] = sort(dataDate);
%updatedDateIdx = find(dataDate > localLastDate);

if ind_dates % individual dates
    dateRange = convertDateStr(dateRange,'file_date');
    X = cellfun(@(c)strcmp(c,dataDateStr),dateRange,'UniformOutput',false);
    date_idx = any(vertcat(X{:}));
    updatedDataFiles = fileList(date_idx);    
else
    if isempty(endDate)
        updatedDataFiles = fileList(sortIdx(dataDate > localLastDate));
    else
        updatedDataFiles = fileList(sortIdx(dataDate >= localLastDate & dataDate <= endDate));
    end
end


% compare existing files
% update if file size and date is bigger
local_fnl = {local_fileList.name};
source_fnl = {source_fileList.name};

for lcIdx = 1:length(local_fnl)
    srcIdx = find(strcmp(local_fnl{lcIdx}, source_fnl));
    if ~isempty(srcIdx)
        if source_fileList(srcIdx).bytes > local_fileList(lcIdx).bytes
            if source_fileList(srcIdx).datenum > local_fileList(lcIdx).datenum
                updatedDataFiles{length(updatedDataFiles)+1} = local_fileList(lcIdx).name;
            end
        end
    end
end


if ~isempty(updatedDataFiles)
    pfileList = getIncFiles(data_source, updatedDataFiles, MSN, sbjIDs);  
    % copy data files to local data folder
    for filenum = 1:length(pfileList)
        sourceData = fullfile(data_source, pfileList{filenum});
        copyfile(sourceData, local_data);
        fprintf('%s copied\n',pfileList{filenum})
    end
else
    fprintf('no updated file\n')
end
end

function mfiles = getIncFiles(data_source, updatedDataFiles, MSN, sbjIDs)

% quickly parse MSN only from file list
msnListByFile = cell(length(updatedDataFiles),35);
sbjListByFile = cell(length(updatedDataFiles),35);
for filenum = 1:length(updatedDataFiles)
    fString = fileread(fullfile(data_source,char(updatedDataFiles(filenum)))); % read the file as one string
    msns = regexp(fString,'MSN:\s([\w-]+)','tokens');
    subjects = regexp(fString,'Subject:\s([\w-\s]*)\r','tokens');
    %msnListByFileNum = 1;
    for msnNum = 1:min(length(msns), length(subjects))
        msnListByFile(filenum,msnNum) = msns{msnNum};
        sbjListByFile(filenum,msnNum) = subjects{msnNum};
        %msnListByFileNum = msnListByFileNum + 1;
    end
end
% handle regexp subject IDs
if ~iscell(sbjIDs) && ischar(sbjIDs) && ~isempty(sbjIDs)
    matchSbj = regexp(sbjListByFile, sbjIDs,'match');
    sbjIDs = unique([matchSbj{:}]);
end

if iscell(sbjIDs)&& ~isempty(sbjIDs)
    sbjCheck = true;
else
    sbjCheck = false;
end
% find file list that have requested MSN & subject IDs
for msnSelected = 1:length(MSN)
    for filenum=1:length(updatedDataFiles)
        tmp_msnbyFile = msnListByFile(filenum,:);
        tmp_msnbyFile = tmp_msnbyFile(~cellfun('isempty',tmp_msnbyFile));

        if sbjCheck
            tmp_sbjbyFile = sbjListByFile(filenum,:);
            tmp_sbjbyFile = tmp_sbjbyFile(~cellfun('isempty',tmp_sbjbyFile));
            sbjIn = false;
            for sbjIdx = 1:length(sbjIDs)
                sbjIn = sbjIn || any(strcmp(sbjIDs(sbjIdx),tmp_sbjbyFile));
            end
        else
            sbjIn = true;
        end
        parsingFiles(msnSelected,filenum) = sbjIn && any(strcmp(MSN{msnSelected},tmp_msnbyFile));
    end
end

mfiles = updatedDataFiles(any(parsingFiles,1));

end