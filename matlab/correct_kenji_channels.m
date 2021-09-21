% correct_kenji_channels
% Kenji channels were flipped. 
% this script is ran following reorder_xml_channels.ipynb
% it writes new mapping to:
%                       basename.session
%                       basename.sessionInfo
%                       basename.cell_metrics
%                       classification_DeepSuperficial

df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');

date_check_1 = '20-Sep-2021';
date_check_2 = '21-Sep-2021';

idx = contains(df.basepath,'Kenji');
basepaths = unique(df.basepath(idx));


for i = 1:length(basepaths)
    
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    disp([num2str(i),'  ',basepath])
    
    % check to see if file was already ran
    if check_date(basepath,basename,date_check_1,date_check_2)
        continue
    end

    main(basepath,basename)
    
    close all
end

function datepass = check_date(basepath,basename,date_check_1,date_check_2)
file = dir(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']));
image_file = dir(fullfile(basepath,'deepSuperficial_classification_fromRipples.png'));

datepass = [contains(image_file.date,date_check_1) || contains(image_file.date,date_check_2)] &&...
        [contains(file.date,date_check_1) || contains(file.date,date_check_2)];
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
        session.extracellular.nElectrodeGroups = length(session.extracellular.electrodeGroups.channels);
        session.extracellular.nSpikeGroups = length(session.extracellular.electrodeGroups.channels);
    end
    % update basename.sessionInfo
    sessionInfo.AnatGrps = sessionInfo_new.AnatGrps;
    sessionInfo.SpkGrps = sessionInfo_new.SpkGrps;
    sessionInfo.ElecGp = {sessionInfo_new.AnatGrps.Channels};
    sessionInfo.nElecGps = session.extracellular.nSpikeGroups;
    
    % update cell_metrics
    cell_metrics.general.electrodeGroups = session.extracellular.electrodeGroups.channels; 
    
    % save modified files
    save(fullfile(basepath,[basename,'.session.mat']),'session')
    save(fullfile(basepath,[basename,'.sessionInfo.mat']),'sessionInfo')
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
    
    % run deep superficial classification 
    classification_DeepSuperficial(session,'basepath',basepath);
end

