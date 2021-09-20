% correct_kenji_channels
% Kenji channels were flipped. 
% this script is ran following reorder_xml_channels.ipynb
% it writes new mapping to:
%                       basename.session
%                       basename.sessionInfo
%                       basename.cell_metrics
%                       classification_DeepSuperficial

df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');

idx = contains(df.basepath,'Kenji');
basepaths = unique(df.basepath(idx));

for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    disp([num2str(i),'  ',basepath])
    
    % check to see if file was already ran
    file = dir(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']));
    if contains(file.date,date)
        continue
    end
    main(basepath,basename)
    
    close all
end
function main(basepath,basename)
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.sessionInfo.mat']))
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))

    sessionInfo_new = LoadXml(fullfile(basepath,[basename, '.xml']));
    
    % update basename.session    
    if isfield(sessionInfo_new,'SpkGrps')
        session.extracellular.spikeGroups.channels = cellfun(@(x) x+1,{sessionInfo_new.SpkGrps.Channels},'un',0); % Spike groups
        session.extracellular.electrodeGroups.channels = cellfun(@(x) x+1,{sessionInfo_new.SpkGrps.Channels},'un',0);
    else
        warning('No spike groups exist in the xml. Anatomical groups used instead')
        session.extracellular.spikeGroups.channels = cellfun(@(x) x+1,{sessionInfo_new.AnatGrps.Channels},'un',0); % Spike groups
        session.extracellular.electrodeGroups.channels = cellfun(@(x) x+1,{sessionInfo_new.AnatGrps.Channels},'un',0);
    end
    % update basename.sessionInfo
    sessionInfo.AnatGrps = sessionInfo_new.AnatGrps;
    sessionInfo.SpkGrps = sessionInfo_new.SpkGrps;
    
    % update cell_metrics
    cell_metrics.general.electrodeGroups = session.extracellular.electrodeGroups.channels; 
    
    % save modified files
    save(fullfile(basepath,[basename,'.session.mat']),'session')
    save(fullfile(basepath,[basename,'.sessionInfo.mat']),'sessionInfo')
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
    
    % run deep superficial classification 
    classification_DeepSuperficial(session,'basepath',basepath);
end

