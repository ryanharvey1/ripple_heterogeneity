
electrodes = 4
waveforms_ = reshape(waveforms, [electrodes,samples/electrodes,length(waveforms)/samples]);

waveforms_ = reshape(waveforms,10,32,[]);

%% custom for A:\Data\GrosmarkAD\Achilles\Achilles_10252013
old_spikes.spikes.region

spikes.region = spikes.shankID;
spikes.region = repmat({'unknown'},length(spikes.shankID),1);

spikes.region(spikes.shankID <= 6) = {'lCA1'}
spikes.region(spikes.shankID > 6) = {'rCA1'}

spikes.region = spikes.region';

save(fullfile(basepath,[basename '.spikes.cellinfo.mat']),'spikes')

%%
load('ripCh.mat')
swrCh.ripple = ripCh;
ripples = bz_DetectSWR([swrCh.ripple swrCh.sharpwave],'saveMat',true);

%%
basepath = 'A:\Data\AB3\AB3_47_49';
basename = bz_BasenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

load('cell_waveforms.mat');
max(spikes.filtWaveform{UID}) - min(spikes.filtWaveform{UID});
cell_metrics.peakVoltage

cell_metrics = ProcessCellMetrics('basepath',basepath,...
    'showGUI',false,...
    'spikes',spikes,...
    'getWaveformsFromDat',false);

%% custom bz_getRipSpikes for A:\OptoMECLEC\OML18 & OML19
% custom because some session contain pre/task/task/post & other combos

basepath = 'A:\OptoMECLEC\OML19\day3';
basename = bz_BasenameFromBasepath(basepath);
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']))

cell_metrics.CA1depth = cell_metrics.deepSuperficialDistance;
save(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']),'cell_metrics')

load(fullfile(basepath,[basename '.ripples.events.mat']))

load(fullfile(basepath,[basename '.session.mat']))

for i= 1:length(session.epochs)
    disp(session.epochs{i}.name)
end

epochs = {'pre','task','post'};
epoch_struct.pre = [session.epochs{1}.startTime,session.epochs{1}.stopTime];
epoch_struct.task = [session.epochs{2}.startTime,session.epochs{2}.stopTime];
epoch_struct.post = [session.epochs{3}.startTime,session.epochs{3}.stopTime];

SWRunitMetrics = [];
SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
                                            ripples,...
                                            SWRunitMetrics,...
                                            epochs,...
                                            epoch_struct,...
                                            'pre',...
                                            'pre'); 
                                        
SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
                                            ripples,...
                                            SWRunitMetrics,...
                                            epochs,...
                                            epoch_struct,...
                                            'task',...
                                            'task'); 
                                        
SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
                                            ripples,...
                                            SWRunitMetrics,...
                                            epochs,...
                                            epoch_struct,...
                                            'post',...
                                            'post');     
                                        
save(fullfile(basepath,[basename '.SWRunitMetrics.mat']),'SWRunitMetrics')
                                        
%%

function SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
                                                    ripples,...
                                                    SWRunitMetrics,...
                                                    epochs,...
                                                    epoch_struct,...
                                                    epoch_search,...
                                                    epoch_label)
                                                
start_end = epoch_struct.(epochs{contains(epochs,epoch_search)});
if isfield(ripples,'times')
    idx = ripples.times(:,1) >= start_end(1) &...
        ripples.times(:,2) <= start_end(2);
    ripSpk = bz_getRipSpikes('basepath',basepath,...
                        'events',ripples.times(idx,:),...
                        'saveMat',false); 
elseif isfield(ripples,'timestamps')
    idx = ripples.timestamps(:,1) >= start_end(1) &...
        ripples.timestamps(:,2) <= start_end(2);
    ripSpk = bz_getRipSpikes('basepath',basepath,...
                        'events',ripples.timestamps(idx,:),...
                        'saveMat',false); 
end

SWRunitMetrics.(epoch_label) = bz_unitSWRmetrics(ripSpk);

end

