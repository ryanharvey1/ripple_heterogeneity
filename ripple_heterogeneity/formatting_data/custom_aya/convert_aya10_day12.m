basepath = 'Z:\Data\AYAold\AYA10\day12'
basename = basenameFromBasepath(pwd);
load('day12.session.mat')
load('cell_waveforms.mat');
load('day12.av_waveform.mat')

[~,midx] = min(squeeze(waveforms(1,:,:)));
timeWaveform = ((1:size(waveforms,2)) - mode(midx)) / session.extracellular.sr * 1000;

for shank_idx = 1:length(av_waveform)
    UID = spikes.UID(spikes.shankID == shank_idx);
    
    for cell_i = 1:length(av_waveform{shank_idx})
        max_ = max(abs(av_waveform{shank_idx}{cell_i}),[],2);
        [~,idx] = max(max_);
        spikes.filtWaveform{spikes.UID == UID(cell_i)} =...
            av_waveform{shank_idx}{cell_i}(idx,:);
        
        spikes.maxWaveformCh1(spikes.UID == UID(cell_i)) =...
            session.extracellular.electrodeGroups.channels{shank_idx}(idx);
        
        spikes.maxWaveformCh(spikes.UID == UID(cell_i)) = spikes.maxWaveformCh1(spikes.UID == UID(cell_i))-1;
        
        spikes.peakVoltage(spikes.UID == UID(cell_i)) =...
            double(range(spikes.filtWaveform{spikes.UID == UID(cell_i)}));
        
        spikes.filtWaveform_all{spikes.UID == UID(cell_i)} = ...
            av_waveform{shank_idx}{cell_i};
        
        spikes.filtWaveform_all_std{spikes.UID == UID(cell_i)} =...
            std(av_waveform{shank_idx}{cell_i});
        
        spikes.channels_all{spikes.UID == UID(cell_i)} = session.extracellular.electrodeGroups.channels{shank_idx};
        
        spikes.timeWaveform{spikes.UID == UID(cell_i)} = timeWaveform;
    end
end
save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')
%%

load('posTrials.mat')

% idx = find(~isnan(posTrials{1, 1}(:,2)))
% posTrials{1, 1}(idx(1),1)
%
% offset = 0;
% tx = []
% for i = 1:length(posTrials)
%     posTrials{i}(:,1)
% end

% 0,1640, tmaze
% 1640,8537, sleep
% 8537,10562, tmaze
% 10562,14307.094, sleep

%%

cell_metrics = ProcessCellMetrics('basepath',basepath,...
    'showGUI',true,...
    'spikes',spikes,...
    'getWaveformsFromDat',false,...
    'manualAdjustMonoSyn',false,...
    'session',session,'forceReload',true);

%%
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "Hoja1";
opts.DataRange = "A1:G77";

% Specify column names and types
opts.VariableNames = ["VarName1", "QUALITY", "REGION", "LAYER", "TYPE", "SHANK", "SUBTYPE"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "string"];

% Specify variable properties
opts = setvaropts(opts, "SUBTYPE", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "SUBTYPE", "EmptyFieldRule", "auto");

% Import the data
whowhat = readtable("Z:\Data\AYAold\AYA10\day12\day12.xlsx", opts, "UseExcel", false);

whowhat(1,:) = [];
% Clear temporary variables
clear opts
%%
whowhat = readwhowhat("Z:\Data\AYAold\AYA10\day12\day12.xlsx");

for i = 1:numel(spikes.ts)
    try
        if whowhat.region(i) == 1
            cell_metrics.brainRegion{i} = 'CA1';
        elseif whowhat.region(i) == 2
            cell_metrics.brainRegion{i} = 'CA2';
        elseif whowhat.region(i) == 3
            cell_metrics.brainRegion{i} = 'CA3';
        elseif whowhat.region(i) == 4
            cell_metrics.brainRegion{i} = 'DG';
        else
            cell_metrics.brainRegion{i} = 'Unknown';
        end
    catch
        cell_metrics.brainRegion{i} = 'LEC';
    end
end

save([basename '.cell_metrics.cellinfo.mat'],'cell_metrics');

%%

%  ripples = FindRipples(pwd,swrCh.Ripple_Channel,'noise',swrCh.Noise_Channel,'saveMat',true);
 ripples = DetectSWR([54 60],'saveMat',true);
