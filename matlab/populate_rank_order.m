% add_rank_order_to_cell_metrics

df = readtable('D:\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);


WaitMessage = parfor_wait(length(basepaths));
parfor i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    main(basepath,basename)
    WaitMessage.Send;
end
WaitMessage.Destroy;

function main(basepath,basename)
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))

if isfield(cell_metrics,'rankorder')
    return
end
load(fullfile(basepath,[basename,'.session.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
parameters = cell_metrics.general.processinginfo.params;

cell_metrics = customCalculations.('rank_order_metrics')(cell_metrics,session,{spikes},parameters);

save_new_metrics(basepath,basename,cell_metrics)

end
function save_new_metrics(basepath,basename,cell_metrics)
dirname = 'revisions_cell_metrics';
saveAs = 'cell_metrics';
if ~(exist(fullfile(basepath,dirname),'dir'))
    mkdir(fullfile(basepath,dirname));
end
if exist(fullfile(basepath,[basename,'.',saveAs,'.cellinfo.mat']),'file')
    copyfile(fullfile(basepath,[basename,'.',saveAs,'.cellinfo.mat']),...
        fullfile(basepath, dirname, [saveAs, '_',datestr(clock,'yyyy-mm-dd_HHMMSS'), '.mat']));
end
save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end
% WaitMessage = parfor_wait(length(basepaths));
% for i = 1:length(basepaths)
%     basepath = basepaths{i};
%     basename = basenameFromBasepath(basepath);
%     
%     disp(basepath)
%     
%     if exist([basepath filesep basename '.rankStats.mat'],'file')
%         continue
%     end
% 
%     load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
%     load(fullfile(basepath,[basename,'.ripples.events.mat']))
%     
%     ripSpk = getRipSpikes('basepath',basepath,...
%                     'spikes',spikes,...
%                     'events',ripples.timestamps,...
%                     'saveMat',false);
%                 
%     [rankStats] = RankOrder('basepath',basepath,...
%                             'spkEventTimes',ripSpk,...
%                             'saveMat',true,...
%                             'numRep',0,...
%                             'doPlot',false);
%     WaitMessage.Send;
% end
% WaitMessage.Destroy; 


% function run_it = check_file(basepath,basename)
% load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
% run_it = ~isfield(cell_metrics,'ripple_particip');
% end




% basepath = 'A:\Data\Kenji\ec013.702_724'
% load(fullfile(basepath,[basenameFromBasepath(basepath),'.spikes.cellinfo.mat']))
% load(fullfile(basepath,[basenameFromBasepath(basepath),'.ripples.events.mat']))
% 
% spkEventTimes = getRipSpikes('basepath',basepath,...
%                     'spikes',spikes,...
%                     'events',ripples.timestamps,...
%                     'saveMat',false);
% 
% [rankStats] = RankOrder('basepath',basepath,'spkEventTimes',ripSpk,'saveMat',true,'numRep',250)
% 
% idx = rankStats.pvalEvents < 0.05;
% nanmedian(rankStats.rankUnits(:,idx)');
% 
% figure;
% scatter(nanmedian(rankStats.rankUnits'),nanmedian(rankStats.rankUnits(:,idx)'))
% xlim([0,1])
% ylim([0,1])
% hold on
% plot([0,1],[0,1])
% 
% 
% rankStats.rankClusters
% 
% 
% get_Rank_matrix
% % Create (#events) x (#units) matrix with position of spikes in event. It
% % considers just the first spikes of each unit
% evtTimes = spkEventTimes.EventRel;
% timeSpike = 'first';
% minUnits = 10;
% normalize = true;
% rankUnits = nan*ones(size(spkEventTimes.UnitEventRel));
% for event = 1:length(evtTimes)
%     % Take into account just first spike
%     if strcmp(timeSpike,'first')
%         units = unique(evtTimes{event}(2,:),'stable');
%     elseif strcmp(timeSpike,'mean')
%         means = [];
%         for jj = unique(evtTimes{event}(2,:))
%             means  = [means [mean(evtTimes{event}(1,evtTimes{event}(2,:)==jj)); jj]];
%         end
%         if ~isempty(means)
%             means = sortrows(means')';
%             units = means(2,:);
%         else
%             units = [];
%         end
%     else
%         warning('The variable "timeSpike" is invalid');
%     end
%     nUnits = length(units);
%     % Set event as nan if it has no enough units
%     if nUnits < minUnits
%         rankUnits(:,event) = nan;
%     % Rank units
%     else
%         rankUnits(units,event) = 1:nUnits;
%         % If normalize, set ranks between 0 and 1        
%         if normalize
%             rankUnits(units,event) = rankUnits(units,event) / nUnits;
%         end
%     end
% end


