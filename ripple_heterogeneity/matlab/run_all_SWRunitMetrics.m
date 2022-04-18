% run_all_SWRunitMetrics
df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepaths = unique(df.basepath);

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.ripples.events.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
    load(fullfile(basepath,[basename,'.SleepState.states.mat']))

    SleepState.ints.NREMstate
    % locate NREMstate epochs (NREMstate is co-localized with theta by def)
    rem_n = find(strcmp(SleepState.idx.statenames,'NREM'));
    idx_rem = ToIntervals(SleepState.idx.states == rem_n);

    % locate awake epochs
    wake_n = find(strcmp(SleepState.idx.statenames,'WAKE'));
    idx_wake = SleepState.idx.states == wake_n;
    
    % locate nontheta epochs
    theta_n = find(strcmp(SleepState.idx.theta_epochs.statenames,'nonTHETA'));
    idx_theta = SleepState.idx.theta_epochs.states == theta_n;

    % locate awake nontheta epochs
    idx_wake_nontheta = ToIntervals(idx_wake & idx_theta);
    
    ripples_awake_nontheta = eventIntervals(ripples,idx_wake_nontheta,true);
    
    ripples_nrem = eventIntervals(ripples,SleepState.ints.NREMstate,true);
    
    ripSpk = bz_getRipSpikes('basepath',basepath,...
                            'basename',basename,...
                            'spikes',spikes,...
                            'events',ripples_nrem.timestamps,...
                            'saveMat',false); 
                        
    avg_rate = spikes.total / (max(cell2mat(spikes.times')) - min(cell2mat(spikes.times')));

    SWRunitMetrics.(epoch_label) = unitSWRmetrics(ripSpk,'baseFR',avg_rate');
                 
%     parse_pre_task_post(session,basepath,basename,ripples,spikes)
    
end

field = fields(ripples);

for i = 1:size(SleepState.ints.NREMstate,1)
    ripples.timestamps > SleepState.ints.NREMstate(i)
    ripples.(field)
end