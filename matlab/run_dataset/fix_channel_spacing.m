%% The issue: all sessions had improper channel spacing,
%   which was propagated into deep/sup pyr dist
%
%   most channels are ~20um apart
%
% get data sheet
if exist('Z:\home\ryanh\projects\ripple_heterogeneity\probe_df.csv','file')
    probe_df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\probe_df.csv');
else
    df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
    df_probe = table();
    for i = 1:length(df.basepath)
        basepath = df.basepath{i};
        basename = basenameFromBasepath(basepath);
        load(fullfile(basepath,[basename,'.session.mat']))
        
        temp_df = table();
        temp_df.basepath = basepath;
        temp_df.animal = animalFromBasepath(basepath);
        temp_df.probesLayout = session.analysisTags.probesLayout;
        try
            temp_df.probesVerticalSpacing = session.analysisTags.probesVerticalSpacing;
        catch
            temp_df.probesVerticalSpacing = NaN;
        end
        
        C{i} = table2cell(temp_df);
    end
    probe_df = cell2table(vertcat(C{:}),...
        "VariableNames",["basepath" "animal" "probesLayout" "probesVerticalSpacing"]);
    writetable(probe_df,'Z:\home\ryanh\projects\ripple_heterogeneity\probe_df.csv')
end

% iter over each session
for i = 1:length(probe_df.basepath)
%     if probe_df.fixed(i) == true
%         continue
%     end
    fix_it(probe_df.basepath{i},probe_df.real_spacing(i))
    probe_df.fixed(i) = true;
    writetable(probe_df,'Z:\home\ryanh\projects\ripple_heterogeneity\probe_df.csv')
end

% func to fix channel dist and prop changes throughout all data types
function fix_it(basepath,real_spacing)

    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']))
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

    for group_i = 1:length(deepSuperficialfromRipple.ripple_channels)
        % locate channels in group
        channels = deepSuperficialfromRipple.ripple_channels{group_i};
        [~,b] = ismember(channels,deepSuperficialfromRipple.channel);
        % find current offset
        % sort to account for mapping differences
        sorted_channels = sort(deepSuperficialfromRipple.channelDistance(b));
        current_channel_offset = mode(diff(sorted_channels));
        % find ratio between real and current
        ratio = real_spacing / current_channel_offset;
        % adjust current by ratio to make real channel dist
        deepSuperficialfromRipple.channelDistance(b) =...
            deepSuperficialfromRipple.channelDistance(b) * ratio;
    end
    
    % assign spacing meta data, this is hidden all over the place
    % in session
    session.animal.probeImplants{1, 1}.verticalSpacing = real_spacing;
    session.animal.probeImplants{1, 1}.supplier = [];
    session.animal.probeImplants{1, 1}.descriptiveName = [];
    session.analysisTags.probesVerticalSpacing = real_spacing;
    % in deepSuperficialfromRipple file
    deepSuperficialfromRipple.processinginfo.params.verticalSpacing = real_spacing;
    % in cell metrics
    cell_metrics.general.SWR = deepSuperficialfromRipple;
    cell_metrics.general.SWR_batch = deepSuperficialfromRipple;
    cell_metrics.general.processinginfo.deepSuperficialfromRipple.params.verticalSpacing = real_spacing;
    cell_metrics.general.deepSuperficial_file =...
        fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']);
    %     cell_metrics.general.processinginfo.params.session.analysisTags.probesVerticalSpacing = real_spacing;
    cell_metrics.general.processinginfo.params.session = session;
    
    % for some reason, the basename can change... fix this here
    if ~contains(cell_metrics.general.basename,basename)
        cell_metrics.general.basename = basename;
    end
    
    % update cell locations
    for j = 1:cell_metrics.general.cellCount
        if isnan(spikes.maxWaveformCh1(j))
            cell_metrics.deepSuperficial(j) = {'Unknown'};
            cell_metrics.deepSuperficialDistance(j) = NaN;
        else
            cell_metrics.deepSuperficial(j) =...
                deepSuperficialfromRipple.channelClass(spikes.maxWaveformCh1(j));
            cell_metrics.deepSuperficialDistance(j) =...
                deepSuperficialfromRipple.channelDistance(spikes.maxWaveformCh1(j));
        end
    end
    
    save(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),'deepSuperficialfromRipple')
    save(fullfile(basepath,[basename,'.session.mat']),'session')
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end

%         [a,b] = ismember(channels,deepSuperficialfromRipple.channel);
%         channel_idx = find(deepSuperficialfromRipple.channelDistance(b) == 0);
%         if ~isempty(channel_idx)
%             new_dist = linspace(0,length(b)*real_spacing,length(b)) - channel_idx*real_spacing;
%             deepSuperficialfromRipple.channelDistance(b) = new_dist;
%         else
%             ratio = real_spacing / mode(diff(deepSuperficialfromRipple.channelDistance(b)))
%             deepSuperficialfromRipple.channelDistance(b) =...
%                 deepSuperficialfromRipple.channelDistance(b) * ratio;
%         end
% ccf.x = zeros(16,1); % Setting x-position to zero
% ccf.y = zeros(16,1); % Setting y-position to zero
% ccf.z = -20*[0:15]'; % Assigning negative depth values, starting at y=0 for channel 1





