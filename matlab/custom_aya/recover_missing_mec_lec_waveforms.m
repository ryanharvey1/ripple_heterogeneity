% recover_missing_mec_lec_waveforms

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepaths = unique(df(contains(df.animal,'AYA'),:).basepath);


for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    try
        load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
%         load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
        load(fullfile(basepath,[basename, '.session.mat']));
    catch
        continue
    end
    if any(isnan(spikes.maxWaveformCh1))
        
        files = dir(fullfile(basepath,'clusteringG','**','*.av_waveform.mat'));
        for sk = 1:length(files)
            load(fullfile(files(sk).folder,files(sk).name))
            shank_idx = str2double(extractAfter(files(sk).folder,'shk'));            
            spikes.filtWaveform_all(spikes.shankID == shank_idx) = av_waveform{shank_idx};
            
            spikes.filtWaveform(spikes.shankID == shank_idx) = av_waveform{shank_idx};
            
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
                
                spikes.filtWaveform_all_std{spikes.UID == UID(cell_i)} =...
                    std(av_waveform{shank_idx}{cell_i});
            end
        end
        disp([files.folder])
        disp(spikes.maxWaveformCh)
        disp(basepath)
        disp('dbcont')
        keyboard
        save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')
        channel_mapping('basepath',basepath)
        

    end
end

for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    try
        load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    catch
        continue
    end
    
    if any(isnan(cat(2,cell_metrics.waveforms.filt{:})))
        load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
        load(fullfile(basepath,[basename, '.session.mat']));

        for s = 1:length(spikes.times)
            spikes.times{s} = unique(spikes.times{s});
            spikes.total(s) = length(spikes.times{s});
        end
        spikes = get_spindices(spikes);
        save(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'spikes')

        cell_metrics = ProcessCellMetrics('basepath',basepath,...
            'showGUI',false,...
            'spikes',spikes,...
            'getWaveformsFromDat',false,...
            'manualAdjustMonoSyn',false,...
            'session',session,'forceReload',true);
        
        save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
    end
%     if any(isnan(cell_metrics.waveforms.filt_zscored))
%         load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
%         load(fullfile(basepath,[basename, '.session.mat']));
% 
%         cell_metrics = ProcessCellMetrics('basepath',basepath,...
%             'showGUI',false,...
%             'spikes',spikes,...
%             'getWaveformsFromDat',false,...
%             'manualAdjustMonoSyn',false,...
%             'session',session,'forceReload',true);
%         
%         save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
%     end
end