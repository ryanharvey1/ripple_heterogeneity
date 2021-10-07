% add_epochs
% see if basename.session or basename.cell_metrics.cellinfo are missing
% epoch data...if so add it from either .behavEpochs.mat or '_sessInfo.mat'
%

% read in dataframe (here, this .csv only needs to have df.basepath)
df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');

% find unique basepaths to iter over
basepaths = unique(df.basepath);

% iter over basepaths
for i = 1:length(basepaths)
    disp(basepaths{i})
    basename = basenameFromBasepath(basepaths{i});
    load(fullfile(basepaths{i},[basename,'.cell_metrics.cellinfo.mat']))
    load(fullfile(basepaths{i},[basename,'.session.mat']))
    
    % check for epochs elsewhere
    if ~isfield(cell_metrics.general,'epochs')
        cell_metrics.general.epochs = session.epochs;
        save(fullfile(basepaths{i},[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
    end
    if length(cell_metrics.general.epochs)==1
        % try to find epochs in behavEpochs (AYA)
        if exist(fullfile(basepaths{i},[basename,'.behavEpochs.mat']),'file')
            load(fullfile(basepaths{i},[basename,'.behavEpochs.mat']))
            
            % get epoch labels
            epochs = fields(behavEpochs.int);
            
            % check consistency of behavEpochs
            behavEpochs = check_and_fix_behavEpochs(epochs,behavEpochs);
            
            % get epochs again in case they changed
            epochs = fields(behavEpochs.int);
            
            % iter through each epoch and add to .session & cell metrics
            for e = 1:length(epochs)
                if isempty(behavEpochs.int.(epochs{e}))
                    continue
                end
                cell_metrics.general.epochs{e}.name = epochs{e};
                cell_metrics.general.epochs{e}.startTime = behavEpochs.int.(epochs{e})(1);
                cell_metrics.general.epochs{e}.stopTime = behavEpochs.int.(epochs{e})(2);
                
                session.epochs{e}.name = epochs{e};
                session.epochs{e}.startTime = behavEpochs.int.(epochs{e})(1);
                session.epochs{e}.stopTime = behavEpochs.int.(epochs{e})(2);
            end
        elseif exist(fullfile(basepaths{i},[basename,'_sessInfo.mat']),'file')
            % try to find epochs in _sessInfo (\GrosmarkAD)
            load(fullfile(basepaths{i},[basename,'_sessInfo.mat']))
            epochs = fields(sessInfo.Epochs);
            epochs = epochs(contains(epochs,'Epoch'));
            for e = 1:length(epochs)
                if isempty(sessInfo.Epochs.(epochs{e}))
                    continue
                end
                cell_metrics.general.epochs{e}.name = epochs{e};
                cell_metrics.general.epochs{e}.startTime = sessInfo.Epochs.(epochs{e})(1);
                cell_metrics.general.epochs{e}.stopTime = sessInfo.Epochs.(epochs{e})(2);
                
                session.epochs{e}.name = epochs{e};
                session.epochs{e}.startTime = sessInfo.Epochs.(epochs{e})(1);
                session.epochs{e}.stopTime = sessInfo.Epochs.(epochs{e})(2);
            end
        end
        save(fullfile(basepaths{i},[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
        save(fullfile(basepaths{i},[basename,'.session.mat']),'session')
    end
end

function behavEpochs = check_and_fix_behavEpochs(epochs,behavEpochs)
for e = 1:length(epochs)
    % find if single epoch has > 1 epoch
    if size(behavEpochs.int.(epochs{e})) ~= [1,1]
        time_stamps = [];
        epoch = [];
        % unpack each epoch
        for e = 1:length(epochs)
            epoch = [epoch;repmat({epochs{e}},size(behavEpochs.int.(epochs{e}),1),1)];
            time_stamps = [time_stamps;behavEpochs.int.(epochs{e})];
        end
        % make sure epochs are in order
        [~,idx] = sort(time_stamps(:,1));
        epoch = epoch(idx);
        time_stamps = time_stamps(idx,:);
        
        % check if epoch labels have unique names
        % if the label comes up again, put the iter number after it
        for e = 1:length(epoch)
            locs = find(contains(epoch,epoch{e}));
            if length(locs) > 1
                for l = 2:length(locs)
                    epoch{locs(l)} = [epoch{locs(l)},'_',num2str(l)];
                end
            end
        end
        % clear and add epochs
        behavEpochs.int = [];
        for e = 1:length(epoch)
            behavEpochs.int.(epoch{e}) = time_stamps(e,:);
        end
        return
    end
end
end
