% update_or_opto_stim

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);
basepaths = basepaths(contains(basepaths,'ORproject'));

for basepath = basepaths'
    basepath = basepath{1};
    disp(basepath)
    main(basepath)
end

function main(basepath)
    basename = basenameFromBasepath(basepath);
    if exist(fullfile(basepath,[basename,'.optoStim.manipulation.mat']),'file')
        return
    end
    
    old_pulses = load(fullfile(basepath,[basename,'.stimPulses.mat']));
    load(fullfile(basepath,[basename,'.session.mat']));

    old_pulses = sort([old_pulses.pulsesAM;old_pulses.pulsesPM]);
    
    % added closed loop and delayed pulses, pulses lasted for 100ms
    pulses.timestamps(:,1) = old_pulses;
    pulses.timestamps(:,2) = pulses.timestamps(:,1) + 0.1;
    
    % add event label (closed loop is 0 and delayed is 1)
    pulses.eventGroupID = ones(length(pulses.timestamps),1);
    
    % make sure ts are sorted
    [~,idx] = sort(pulses.timestamps(:,1));
    pulses.timestamps = pulses.timestamps(idx,:);
    
    pulses.amplitude = nan(length(pulses.timestamps),1);
    
    pulses.duration = pulses.timestamps(:,2) - pulses.timestamps(:,1);
    
    optoStim.timestamps = pulses.timestamps;
    optoStim.peaks = median(pulses.timestamps,2);
    optoStim.amplitude = pulses.amplitude;
    optoStim.amplitudeUnits = 'au';
    optoStim.eventID = pulses.eventGroupID;
    optoStim.eventIDlabels = {'closed_loop'};
    optoStim.eventIDbinary = false;
    optoStim.center = median(pulses.timestamps,2);
    optoStim.duration = pulses.duration;
    optoStim.detectorinfo = 'getAnalogPulses_update_or_opto_stim';
    
    saveStruct(optoStim,'manipulation','session',session);
end
