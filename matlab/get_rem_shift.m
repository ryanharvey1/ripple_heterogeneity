function rem_shift_data = get_rem_shift(varargin) 
% get_rem_shift: compares phase locking in awake vs. rem to locate 
% rem shifting pyr units in deep ca1. 
%
% Based on Mizuseki, et al 2011.
% Neurons with <120° or >300° preferred theta phases during REM were
% designated as REM-shifting cells, whereas those between 120° to 300° 
% were designated as nonshifting cells. 
% 
% Ryan H 2021

p = inputParser;
addParameter(p,'basepath',pwd) % path to folder
addParameter(p,'fig',false) % simple debugging/summary figs
addParameter(p,'passband',[6,12]) % theta band
addParameter(p,'lfp',[]) % structure from getLFP with single theta channel
addParameter(p,'spikes',[]) % number of bins in phase hist
addParameter(p,'savemat',true) % save output to basepath
addParameter(p,'numBins',18) % number of bins in phase hist

parse(p,varargin{:})
basepath = p.Results.basepath; 
fig = p.Results.fig; 
passband = p.Results.passband; 
lfp = p.Results.lfp; 
spikes = p.Results.spikes; 
savemat = p.Results.savemat; 
numBins = p.Results.numBins; 

basename = basenameFromBasepath(basepath);

% locates a deep ca1 channel that maximizes theta power
if isempty(lfp)
    lfp = get_deep_ca1_lfp(basepath,passband);
end

% load your spikes
if isempty(spikes)
    spikes = loadSpikes('basepath',basepath,'basename',basename);
end

try
    load(fullfile(basepath,[basename,'.SleepState.states.mat']));
catch
    warning('you must get SleepState...run SleepScoreMaster.m')
end

% locate awake epochs
wake_n = find(strcmp(SleepState.idx.statenames,'WAKE'));
idx_wake = SleepState.idx.states == wake_n;

% locate rem epochs (rem is co-localized with theta by def)
rem_n = find(strcmp(SleepState.idx.statenames,'REM'));
idx_rem = SleepState.idx.states == rem_n;

% locate theta epochs
theta_n = find(strcmp(SleepState.idx.theta_epochs.statenames,'THETA'));
idx_theta = SleepState.idx.theta_epochs.states == theta_n;
  
% locate awake theta epochs
idx_wake_theta = idx_wake & idx_theta;

% get your phase histograms for awake and rem
PhaseLockingData_wake = phaseModulation(spikes,...
                                    lfp,...
                                    passband,...
                                    'saveMat', false,...
                                    'plotting',false,...
                                    'numBins',numBins,...
                                    'basepath',basepath,...
                                    'intervals',ToIntervals(idx_wake_theta),...
                                    'powerThresh',0);
                                
PhaseLockingData_rem = phaseModulation(spikes,...
                                    lfp,...
                                    passband,...
                                    'saveMat', false,...
                                    'plotting',false,...
                                    'numBins',numBins,...
                                    'basepath',basepath,...
                                    'intervals',ToIntervals(idx_rem),...
                                    'powerThresh',0);
  
% count the number of spikes per cell that were used
spk_count_rem_idx = cellfun('length',PhaseLockingData_rem.spkphases) > 50;

% mean phase angle per unit... in degrees
angles = rad2deg(PhaseLockingData_rem.phasestats.m);

% rem shift class
rem_shift = (angles<120 | angles>300) &...
    PhaseLockingData_rem.phasestats.p<0.01 &...
    spk_count_rem_idx;

% non rem shift class
non_rem_shift = (angles>120 & angles<300) &...
    PhaseLockingData_rem.phasestats.p<0.01 &...
    spk_count_rem_idx;

% store results in rem_shift_data

% get the circular distance between awake and rem
rem_shift_data.circ_dist = circ_dist(PhaseLockingData_rem.phasestats.m,...
                                PhaseLockingData_wake.phasestats.m);                        
rem_shift_data.rem_shift = rem_shift;
rem_shift_data.non_rem_shift = non_rem_shift;
rem_shift_data.PhaseLockingData_rem = PhaseLockingData_rem;
rem_shift_data.PhaseLockingData_wake = PhaseLockingData_wake;
rem_shift_data.detectorParams = p.Results;
rem_shift_data.detectorParams.channels = lfp.channels;
rem_shift_data.detectorParams.samplingRate = lfp.samplingRate;

% 
if savemat
    save(fullfile(basepath,[basename,'.theta_rem_shift.mat']),'rem_shift_data')
end

if fig
    figure;
    states = fields(SleepState.ints);
    viewwin = [SleepState.detectorinfo.StatePlotMaterials.t_clus(1),...
        SleepState.detectorinfo.StatePlotMaterials.t_clus(end)];
    colors = lines(length(states));
    for i = 1:length(states)
        disp(states{i})
        plot(SleepState.ints.(states{i})',(-i)*ones(size(SleepState.ints.(states{i})))',...
            'color',colors(i,:),'LineWidth',8)
        hold on;
    end
    xlim([viewwin]) 
    ylim([-i-1 0])
    set(gca,'YTick',[-i:-1])
    set(gca,'YTickLabel',flipud(states))

    figure;
    x = rad2deg([PhaseLockingData_rem.phasebins;...
        PhaseLockingData_rem.phasebins+2*pi]);
    y = [PhaseLockingData_rem.phasedistros(:,rem_shift);...
        PhaseLockingData_rem.phasedistros(:,rem_shift)];
    plot(x,y,'r')
    hold on
    y = [PhaseLockingData_rem.phasedistros(:,non_rem_shift);...
        PhaseLockingData_rem.phasedistros(:,non_rem_shift)];
    plot(x,y,'k')


    figure;
    x = rad2deg([PhaseLockingData_rem.phasebins;...
        PhaseLockingData_rem.phasebins+2*pi]);
    y = [PhaseLockingData_rem.phasedistros(:,rem_shift);...
        PhaseLockingData_rem.phasedistros(:,rem_shift)];
    plot(x,y,'r')
    hold on
    y = [PhaseLockingData_wake.phasedistros(:,rem_shift);...
        PhaseLockingData_wake.phasedistros(:,rem_shift)];
    plot(x,y,'k')

    figure;
    x = rad2deg([PhaseLockingData_rem.phasebins;...
        PhaseLockingData_rem.phasebins+2*pi]);
    y = [PhaseLockingData_rem.phasedistros(:,rem_shift);...
        PhaseLockingData_rem.phasedistros(:,rem_shift)];
    subplot(2,1,1)
    imagesc(y')
    title('rem')

    hold on
    y = [PhaseLockingData_wake.phasedistros(:,rem_shift);...
        PhaseLockingData_wake.phasedistros(:,rem_shift)];
    subplot(2,1,2)
    imagesc(y')
    title('awake')
end       
end

function lfp = get_deep_ca1_lfp(basepath,passband)
% get_deep_ca1_lfp: locates a deep ca1 channel that maximizes theta power

basename = basenameFromBasepath(basepath);

% load cell_metrics to locate deep ca1 channels
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))

% find deep ca1 channels to check
ca1_idx = contains(lower(cell_metrics.brainRegion),'ca1');
deep_idx = cell_metrics.CA1depth>0;
deep_channels = unique(cell_metrics.maxWaveformCh1(deep_idx & ca1_idx));

if isempty(deep_channels)
    disp('deep channel not found...')
    disp('manually assign deep_channels variable and enter "dbcont"')
    keyboard
end
% load deep channels
lfp = getLFP(deep_channels,'basepath',basepath,'basename',basename);

% get theta power to choose channel
try
    pBand = bandpower(double(lfp.data),...
                    lfp.samplingRate,passband);
                
    pTot = bandpower(double(lfp.data),...
                    lfp.samplingRate,...
                    [1,(lfp.samplingRate/2)-1]);
catch
    for c = 1:size(lfp.data,2)
        pBand(c) = bandpower(double(lfp.data(:,c)),...
                        lfp.samplingRate,passband);
                    
        pTot(c) = bandpower(double(lfp.data(:,c)),...
                            lfp.samplingRate,...
                            [1,(lfp.samplingRate/2)-1]);
    end
end
% find max theta power, normalized by wide band
[~,c_idx] = max(pBand./pTot);

% only leave theta channel
lfp.data = lfp.data(:,c_idx);
lfp.channels = lfp.channels(:,c_idx);
end