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
    % pull in old data container
    old_pulses = load(fullfile(basepath,[basename,'.stimPulses.mat']));
    % load associated behavior stuct which has metadata for stimPulses
    load(fullfile(basepath,[basename,'.behavior.mat']));
    % withing behavior and stimPulses, there are AM and PM pulses
    % if protocol is 1, closed loop ripple stim occured, if 2, then delayed
    % stim occured
    if isfield(old_pulses,'pulsesAM') && isfield(old_pulses,'pulsesPM')
        old_pulses.pulses = sort([old_pulses.pulsesAM;old_pulses.pulsesPM]);
    end
    % some sessions are non standard like:
    %   Z:\Data\ORproject\OR15\hc280118, Z:\Data\ORproject\OR15\hc300118
    if ~isfield(behavior,'sessionAM')
        return
    end
    pulsesD = [];
    pulsesCL = [];
    if behavior.sessionAM.protocol == 1
        pulsesCL = Restrict(old_pulses.pulses,behavior.sessionAM.int);
    elseif behavior.sessionAM.protocol == 2
        pulsesD = Restrict(old_pulses.pulses,behavior.sessionAM.int);
    end
    if behavior.sessionPM.protocol == 1
        pulsesCL = Restrict(old_pulses.pulses,behavior.sessionPM.int);
    elseif behavior.sessionPM.protocol == 2
        pulsesD = Restrict(old_pulses.pulses,behavior.sessionPM.int);
    end
    % put pulses together
    old_pulses = [pulsesCL;pulsesD];
    
    % added closed loop and delayed pulses, pulses lasted for 100ms
    pulses.timestamps(:,1) = old_pulses;
    pulses.timestamps(:,2) = pulses.timestamps(:,1) + 0.1;
    
    % add event label (closed loop is 0 and delayed is 1)
    pulses.eventGroupID = [zeros(length(pulsesCL),1);ones(length(pulsesD),1)];
    
    % not necessary, but make sure ts are sorted
    [~,idx] = sort(pulses.timestamps(:,1));
    pulses.timestamps = pulses.timestamps(idx,:);
    pulses.eventGroupID = pulses.eventGroupID(idx,:);
    
    pulses.amplitude = nan(length(pulses.timestamps),1);
    
    pulses.duration = pulses.timestamps(:,2) - pulses.timestamps(:,1);
    
    optoStim.timestamps = pulses.timestamps;
    optoStim.peaks = median(pulses.timestamps,2);
    optoStim.amplitude = pulses.amplitude;
    optoStim.amplitudeUnits = 'au';
    optoStim.eventID = pulses.eventGroupID;
    optoStim.eventIDlabels = {'closed_loop','delayed'};
    optoStim.eventIDbinary = true;
    optoStim.center = median(pulses.timestamps,2);
    optoStim.duration = pulses.duration;
    optoStim.detectorinfo = 'getAnalogPulses_update_or_opto_stim';
    
    % save stuct as manipulation file
    load(fullfile(basepath,[basename,'.session.mat']));
    saveStruct(optoStim,'manipulation','session',session);
end
