% add_ripples

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\mouse_sessions.csv');
for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    ripple_file = fullfile(basepath,[basename,'.ripples.events.mat']);
    if exist(ripple_file,'file')
        continue
    end
    spikes = importSpikes('basepath',basepath,'cellType', "Pyramidal Cell", 'brainRegion', "CA1");
    if isempty(spikes.times) || length(spikes.UID) < 5
        continue
    end
    load(fullfile(basepath,[basename,'.session.mat']));

    lfp = getLFP(session.brainRegions.CA1.channels,'basepath',basepath,...
        'basename',basename,'noPrompts',true);

    try
        pBand = bandpower(double(lfp.data),lfp.samplingRate,[100,250]);
        pTot = bandpower(double(lfp.data),lfp.samplingRate,[1,(lfp.samplingRate/2)-1]);
    catch
        for c = 1:size(lfp.data,2)
            pBand(c) = bandpower(double(lfp.data(:,c)),lfp.samplingRate,[100,250]);
            pTot(c) = bandpower(double(lfp.data(:,c)),lfp.samplingRate,[1,(lfp.samplingRate/2)-1]);
        end
    end
    ripple_pow = pBand./pTot;
    [~,c_idx] = max(ripple_pow);

    ripple_ch = session.brainRegions.CA1.channels(c_idx);
    ripples = FindRipples(basepath,ripple_ch,'saveMat',true);
    
    ripples = DetectSWR([31 3],'basepath',basepath,'saveMat',true,...
        'thresSDswD',[0.5 2.5],'thresSDrip', [0.25 1],...
        'forceDetect',true,'check',true);

    ripples = eventSpikingTreshold(ripples,'basepath',basepath,'spikes',spikes,'spikingThreshold',0.0001);
    
    save(fullfile(basepath,[basename,'.ripples.events.mat']),'ripples');

end