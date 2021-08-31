% custom_for_girardeau

% has cell metrics but no .spk.
% files = {'A:\Data\GirardeauG\Rat08\Rat08-20130708',...
%     'A:\Data\GirardeauG\Rat08\Rat08-20130709'}


data_path = 'A:\Data\GirardeauG';
tic
files = dir([data_path,'\**\*.spikes.cellinfo.mat']);
disp(toc)


for i = 1:length(files)
    
    basepath = fullfile(files(i).folder);
    disp(basepath)
    
    % Rat07 may not be good. pass it here
    if contains(basepath,'Rat07')
        continue
    end
    
    basename = basenameFromBasepath(basepath);
    
    % check if spikes file exist, make it if not
    if ~exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'file')
        spikes = loadSpikes(...
            'basepath',basepath,...
            'basename',basename,...
            'saveMat',true,...
            'format','Klustakwik',...
            'getWaveformsFromDat',false,...
            'getWaveformsFromSource',true,...
            'forceReload',false);
    end
    
    % load inportant files here
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
    load(fullfile(basepath,[basename,'.ripples.events.mat']))
    
    if ~exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file')
        % load spikes to make sure they are up-to-date
        spikes = loadSpikes(...
            'basepath',basepath,...
            'basename',basename,...
            'saveMat',true,...
            'format','Klustakwik',...
            'getWaveformsFromDat',false,...
            'getWaveformsFromSource',true,...
            'forceReload',false);
        % make cell metrics
        cell_metrics = ProcessCellMetrics('basepath',basepath,...
            'showGUI',false,...
            'spikes',spikes,...
            'getWaveformsFromDat',false,...
            'manualAdjustMonoSyn',false,...
            'session',session);
    else
        load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    end
    % check if epochs are in cell metrics... if not add them
    % first check if the number of epochs are the same
    remake_session = false;
    if exist(fullfile(basepath,[basename,'.cat.evt']),'file')
        events = LoadEvents(fullfile(basepath,[basename,'.cat.evt']));
        if length(session.epochs) ~= length(events.description)/2
            remake_session = true;
        end
    end
    if ~isfield(cell_metrics.general,'epochs') || remake_session
        events = LoadEvents(fullfile(basepath,[basename,'.cat.evt']));
        
        % try to link existing epoch info with epoch data from cat.evt
        if ~remake_session
            try
                for e = 1:length(session.epochs)
                    times = events.time(contains(events.description,session.epochs{e}.name));
                    session.epochs{e}.startTime = times(1);
                    session.epochs{e}.stopTime = times(2);
                end
                % if the above fails, there may be something wrong with .session
            catch
                % rename existing .session with old
                movefile(fullfile(basepath,[basename,'.session.mat']),fullfile(basepath,[basename,'.session_old.mat']),'f');
                % make new .session
                session = sessionTemplate(basepath,'basename',basename,'showGUI',false);
                % iter through .cat.evt and use those epochs to fill .session
                for e = 1:(length(events.description)/2)
                    times = events.time(e:e+1);
                    session.epochs{e}.startTime = times(1);
                    session.epochs{e}.stopTime = times(2);
                    name = strsplit(events.description{e},'-');
                    session.epochs{e}.name = name{end};
                end
            end
        else
            % rename existing .session with old
            movefile(fullfile(basepath,[basename,'.session.mat']),fullfile(basepath,[basename,'.session_old.mat']),'f');
            % make new .session
            session = sessionTemplate(basepath,'basename',basename,'showGUI',false);
            % iter through .cat.evt and use those epochs to fill .session
            for e = 1:(length(events.description)/2)
                times = events.time(e:e+1);
                session.epochs{e}.startTime = times(1);
                session.epochs{e}.stopTime = times(2);
                name = strsplit(events.description{e},'-');
                session.epochs{e}.name = name{end};
            end
        end
        % save .session
        save(fullfile(basepath,[basename,'.session.mat']),'session')
        
        % run ProcessCellMetrics so new epochs will be added to cell_metrics
        cell_metrics = ProcessCellMetrics('basepath',basepath,...
            'showGUI',false,...
            'spikes',spikes,...
            'getWaveformsFromDat',false,...
            'manualAdjustMonoSyn',false,...
            'session',session);
    end
    % fill brain region field, sometimes this info is in region and not brainRegion
    if isempty(cell_metrics.brainRegion)
        cell_metrics.brainRegion = cell_metrics.region;
        save(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics')
    end
    
    % make SWRunitMetrics by parsing epoch data
    if ~exist(fullfile(basepath,[basename '.SWRunitMetrics.mat']),'file')
        parse_pre_task_post(session,basepath,basename,ripples,spikes)
    end
end

% tic
% files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat*']);
% disp(toc)
%
% data_path = 'A:\Data\GirardeauG';
% tic
% files = dir([data_path,'\**\*.spk.*']);
% disp(toc)
%
% sessions = unique({files.folder});
%
% for i = 1:length(sessions)
%     basepath = sessions{i};
%     basename = bz_BasenameFromBasepath(basepath);
%     disp(basepath)
%
%     % check for needed files
%     if exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'file') &&...
%             exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') &&...
%             exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file') &&...
%             exist(fullfile(basepath,[basename,'.SWRunitMetrics.mat']),'file')
%         continue
%     end
%
%     % make sure there are the same n of res and spk
%     if length(dir([basepath,filesep,'*.res.*'])) ~=...
%             length(dir([basepath,filesep,'*.spk.*']))
%         continue
%     end
%
%     % check and make session.mat
%     if ~exist([basename '.session.mat'],'file')
%         session = sessionTemplate(basepath,'showGUI',false);
%         session.epochs = get_kenji_epochs('basepath',basepath,'basename',basename);
%         save(fullfile(basepath,[basename '.session.mat']),'session')
%     else
%         load(fullfile(basepath,[basename '.session.mat']))
%     end
%
% end

% basepath = 'A:\Data\GirardeauG\Rat09\Rat09-20140324'
%
% spikes = loadSpikes('basepath',basepath);
% load('A:\Data\GirardeauG\Rat09\Rat09-20140324\Rat09-20140324.cell_metrics.cellinfo.mat')
%
%
%
% ripSpk = bz_getRipSpikes('basepath',basepath,...
%                             'basename',basename,...
%                             'spikes',spikes,...
%                             'events',ripples.timestamps(idx,:),...
%                             'saveMat',false);
%

