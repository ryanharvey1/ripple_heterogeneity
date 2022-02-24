
sessions = {'AO10/day14', 'AO10/day15', 'AO10/day16', 'AO10/day20', 'AO10/day23', 'AO10/day26', 'AO10/day27',...
    'AO11/day15', 'AO11/day17', 'AO11/day19', 'AO11/day20', 'AO11/day23', 'AO11/day26', 'AO11/day27', 'AO11/day34',...
    'AO12/day8', 'AO12/day9', 'AO12/day11', 'AO12/day13', 'AO12/day15', 'AO12/day17', 'AO12/day18',...
    'AO13/day10', 'AO13/day11', 'AO13/day13', 'AO13/day14', 'AO13/day15', 'AO13/day17', 'AO13/day19',...
    'AO14/day10', 'AO14/day11', 'AO14/day12', 'AO14/day13', 'AO14/day14', 'AO14/day15',...
    'AO15/day9', 'AO15/day10', 'AO15/day11', 'AO15/day12', 'AO15/day13', 'AO15/day14',...
    'AO16/day11', 'AO16/day12',...
    'AO17/day7', 'AO17/day9', 'AO17/day12',...
    'AO19/day9',...
    'AO20/day11', 'AO20/day13', 'AO20/day15', 'AO20/day18', 'AO20/day21', 'AO20/day25', 'AO20/day28',...
    'AO22/day11', 'AO22/day12', 'AO22/day13', 'AO22/day16',...
    'AO23/day13', 'AO23/day18', 'AO23/day21', 'AO23/day22', 'AO23/day25',...
    'AO24/day10', 'AO24/day11',...
    'AO25/day10', 'AO25/day11', 'AO25/day12', 'AO25/day17',...
    'AO26/day10', 'AO26/day11', 'AO26/day13',...
    'AO27/day9', ...
    'AO28/day11', 'AO28/day14',...
    'AO29/day8', 'AO29/day9',...
    'AO31/day10',...
    'AO33/day8', 'AO33/day11',...
    'AO39/day7', 'AO39/day9', 'AO39/day16', 'AO39/day17'}
% basepaths = fullfile('Y:\SMproject',sessions)';
% 
% clear basenames
% for i = 1:length(basepaths)
%     if exist(basepaths{i},'file')
%         basenames{i,1} = basenameFromBasepath(basepaths{i});
%     else
%         disp(basepaths{i})
%     end
% end
% df = table();
% df.basepath = basepaths;
% df.basename = basenames;
% 
% writetable(df,'Z:\home\ryanh\projects\ripple_heterogeneity\mouse_sessions.csv')


% clear basepaths
clear basenames
for i = 1:length(basepaths)
    basename = basenameFromBasepath(basepaths{i});
    basepath = basepaths{i};
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    if ~isfield(cell_metrics,'brainRegion')
        disp(basepath)
        continue
    end
    
    if length(cell_metrics.brainRegion) ~= length(cell_metrics.UID)
        disp(basepath)
        channel_mapping('basepath',basepath)
%         continue
    end
    
    if length(fields(cell_metrics)) == 1
        disp(basepath)
        continue
    end
    basenames{i} = basename;
end

basepaths = basepaths(~cellfun(@isempty,basenames));
basenames = basenames(~cellfun(@isempty,basenames));
% add_files
% files = dir(['Y:\SMproject','\**\*.cell_metrics.cellinfo.mat']);

% Y:\SMproject\AO29\day13

% pull out basepaths and basenames
% clear basepaths
% clear basenames
% for i = 1:length(files)
%     % check if basename.session is there
%     basepath = files(i).folder;
%     basename = basenameFromBasepath(files(i).folder);
%     if ~exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file')
%         disp(['no cell_metrics ',basepath])
%         continue
%     end
%     
%     load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
%     %     try
%     %         [gdata,~,gidx] = unique(cell_metrics.brainRegion);
%     %     catch
%     %        continue
%     %     end
%     if ~isfield(cell_metrics,'brainRegion')
%         disp(['no brain region ',basepath])
%         continue
%     end
%     if any(cellfun(@isempty,cell_metrics.brainRegion))
%         disp(['empty brain region ',basepath])
%         channel_mapping('basepath',basepath)
%         continue
%     end
%     if exist(fullfile(basepath,[basename,'.session.mat']),'file') && isfield(cell_metrics,'general')
%         basepaths{i} = files(i).folder;
%         basenames{i} = basenameFromBasepath(files(i).folder);
%     end
% end
% basepaths = basepaths(~cellfun(@isempty,basepaths));
% basenames = basenames(~cellfun(@isempty,basenames));
% 
% cell_metrics.sessionName(find(cellfun(@isempty,cell_metrics.brainRegion)))
% 
% basepaths(contains(basenames,'day12'))'

% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

% pull up gui to inspect all units in your project
cell_metrics = CellExplorer('metrics',cell_metrics);

