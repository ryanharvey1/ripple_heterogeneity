function epochs = get_kenji_epochs(varargin)
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'basename',[],@isstr); 
addParameter(p,'dirData','A:\Data\Kenji\',@isstr); 

% Parsing inputs
parse(p,varargin{:})
basepath = p.Results.basepath;
basename = p.Results.basename;
dirData = p.Results.dirData;

if isempty(basename)
    parts = strsplit(basepath,filesep);
    basename = parts{end};
end

dirfile = [dirData, basename, filesep];

% load table with kenji session data
load([dirData '\KenjiData3.mat']);

% pull out this current session
info = Beh(ismember(Beh(:,2),basename),:);

% extract timestamps
ts = cumsum(str2double(info(:,7)));

% assign epochs
for i = 1:size(info,1)
    if i==1
        epochs{i}.name = [info{i,4},'_',info{i,5}];
        epochs{i}.startTime = 0;
        epochs{i}.stopTime = ts(i);
    else
        epochs{i}.name = [info{i,4},'_',info{i,5}];
        epochs{i}.startTime = ts(i-1);
        epochs{i}.stopTime = ts(i);
    end
end

end