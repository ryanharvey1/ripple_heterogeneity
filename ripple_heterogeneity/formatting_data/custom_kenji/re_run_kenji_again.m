% re_run_kenji_again
re_run_if_before = '5-Jan-2022 10:00:00';

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);
basepaths = basepaths(contains(basepaths,'Kenji'));

WaitMessage = parfor_wait(length(basepaths),'Waitbar',false);
% loop through each session
parfor i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    if ~exist(fullfile(basepath,[basename,'_old.xml']),'file')
        WaitMessage.Send;
        continue
    end
    pass_metrics = check_file_date(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),...
        re_run_if_before);
    pass_session = check_file_date(fullfile(basepath,[basename,'.session.mat']),...
        re_run_if_before);
    if pass_metrics == 0 && pass_session == 0
        WaitMessage.Send;
        continue
    end
    run_all(basepath,basename)
    WaitMessage.Send;
end
WaitMessage.Destroy;

function run_all(basepath,basename)

mkdir(fullfile(basepath,'old_files'))
files_to_move = {'.sessionInfo.mat','.spikes.cellinfo.mat',...
    '.cell_metrics.cellinfo.mat','.chanCoords.channelInfo.mat',...
    '.session.mat'};
for file = files_to_move
    try
        movefile(fullfile(basepath,[basename,file{1}]),fullfile(basepath,'old_files'))
    catch
        disp([file,' aleady moved'])
    end
end
% check and make session.mat
session = sessionTemplate(basepath,'showGUI',false);

% fill out bad channels
channels = 1:session.extracellular.nChannels;
idx = ismember(channels,cat(2,session.extracellular.electrodeGroups.channels{:}));
session.channelTags.Bad.channels = channels(~idx);

old_session = load(fullfile(basepath,'old_files',[basename,'.session.mat']));

session.epochs = old_session.session.epochs;

save(fullfile(basepath,[basename '.session.mat']),'session')

% load and process spikes
spikes = loadSpikes(...
    'basepath',basepath,...
    'basename',basename,...
    'saveMat',true,...
    'format','Klustakwik',...
    'getWaveformsFromDat',false,...
    'getWaveformsFromSource',true,...
    'forceReload',false);
% 
if mode(diff(session.extracellular.electrodeGroups.channels{1})) == -1
    for UID = spikes.UID
        idx = size(spikes.filtWaveform_all{UID},1):-1:1;
        
        spikes.filtWaveform_all{UID} = spikes.filtWaveform_all{UID}(idx,:);
        spikes.filtWaveform_all_std{UID} = spikes.filtWaveform_all_std{UID}(idx,:);
        
        shank = spikes.shankID(UID);
        [~,index1] = max(max(spikes.filtWaveform_all{UID}') - min(spikes.filtWaveform_all{UID}'));
        spikes.maxWaveformCh(UID) = session.extracellular.electrodeGroups.channels{shank}(index1)-1; % index 0;
        spikes.maxWaveformCh1(UID) = session.extracellular.electrodeGroups.channels{shank}(index1); % index 1;
        spikes.filtWaveform{UID} = spikes.filtWaveform_all{UID}(index1,:);
        spikes.peakVoltage(UID) = max(spikes.filtWaveform{UID}) - min(spikes.filtWaveform{UID});
        spikes.channels_all{UID} = session.extracellular.electrodeGroups.channels{shank};
    end
    save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')
end

% classification_DeepSuperficial(session,'basepath',basepath);

cell_metrics = ProcessCellMetrics('basepath',basepath,...
    'showGUI',false,...
    'spikes',spikes,...
    'getWaveformsFromDat',false,...
    'manualAdjustMonoSyn',false,...
    'session',session);

old_cell_metrics = load(fullfile(basepath,'old_files',...
    [basename,'.cell_metrics.cellinfo.mat']));

if isfield(old_cell_metrics.cell_metrics,'tags')
    cell_metrics.tags = old_cell_metrics.cell_metrics.tags;
end

% save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')
save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')

channel_mapping()
close all
end
function pass = check_file_date(file,date_)
% check if file was ran before or on date
% true if file creation date is less than date_ or doesn't exist

files = dir(file);
if isempty(files)
   pass = 1; 
   return
end
pass = files.datenum < datenum(date_);

% files = dir(file);
% [Y1, M1, D1, H1, ~, ~] = datevec(files.datenum);
% [Y2, M2, D2, H2, ~, ~] = datevec(datenum(date_));
% if D1==D2
%     pass = H1 < H2;
% else
%     pass = D1 < D2;
% end
end
% channels = cat(2,cell_metrics.general.SWR.ripple_channels{:});
% for j = 1:length(cell_metrics.UID)
%     cell_metrics.deepSuperficial(j) = cell_metrics.general.SWR.channelClass(channels == cell_metrics.maxWaveformCh1(j));
%     cell_metrics.deepSuperficialDistance(j) = cell_metrics.general.SWR.channelDistance(channels == cell_metrics.maxWaveformCh1(j));
% end
% 
% cell_metrics.deepSuperficial(j) = deepSuperficial_ChClass(spikes{spkExclu}.maxWaveformCh1(j)); % cell_deep_superficial OK
% cell_metrics.deepSuperficialDistance(j) = deepSuperficial_ChDistance(spikes{spkExclu}.maxWaveformCh1(j)); % cell_deep_superficial_distance