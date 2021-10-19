function swr_seq_pop_isi(varargin)

p = inputParser;
addParameter(p,'basepath',[]) % single basepath to run 1 session
addParameter(p,'df',[]) % data frame with df.basepath
addParameter(p,'binary_class_variable','deepSuperficial',@ischar) % variable in cell_metrics that has groups
addParameter(p,'grouping_names',{'Deep','Superficial'},@iscell) % group names associated with the above
addParameter(p,'restrict_to_brainregion','CA1',@ischar) % brain region to run on (empty to run all)
addParameter(p,'restrict_to_celltype','Pyramidal Cell',@ischar) % cell class to run on (empty to run all)
addParameter(p,'force_run',false,@islogical) % to overwrite results
addParameter(p,'savepath',[]) % path to save results
addParameter(p,'parallel',true,@islogical) % run over sessions in parallel
addParameter(p,'shuffle',false,@islogical) % jitter spike times to make null dist
addParameter(p,'unique_unit_num',3,@isint) % min number of unique units per group per ripple
addParameter(p,'ripple_duration_restrict',[0.05,inf],@isint) % min number of unique units per group per ripple

parse(p,varargin{:})
basepaths = p.Results.basepath;
df = p.Results.df;
binary_class_variable = p.Results.binary_class_variable;
grouping_names = p.Results.grouping_names;
restrict_to_brainregion = p.Results.restrict_to_brainregion;
restrict_to_celltype = p.Results.restrict_to_celltype;
force_run = p.Results.force_run;
savepath = p.Results.savepath;
parallel = p.Results.parallel;
shuffle = p.Results.shuffle;
unique_unit_num = p.Results.unique_unit_num;
ripple_duration_restrict = p.Results.ripple_duration_restrict;

if isempty(basepaths)
    df = readtable('D:\projects\ripple_heterogeneity\sessions.csv');
    basepaths = unique(df.basepath);
end

WaitMessage = parfor_wait(length(basepaths));
if parallel
    parfor i = 1:length(basepaths)
        main(basepaths{i},binary_class_variable,grouping_names,...
            restrict_to_brainregion,restrict_to_celltype,savepath,...
            force_run,shuffle,unique_unit_num,ripple_duration_restrict)
        WaitMessage.Send;
    end
else
    for i = 1:length(basepaths)
        main(basepaths{i},binary_class_variable,grouping_names,...
            restrict_to_brainregion,restrict_to_celltype,savepath,...
            force_run,shuffle,unique_unit_num,ripple_duration_restrict)
        WaitMessage.Send;
    end
end
WaitMessage.Destroy;
end

function main(basepath,binary_class_variable,grouping_names,...
    restrict_to_brainregion,restrict_to_celltype,savepath,force_run,...
    shuffle,unique_unit_num,ripple_duration_restrict)

disp(basepath)

file_str = strsplit(basepath,filesep);
savepath = fullfile(savepath,[file_str{end-1},'_',file_str{end},'.csv']);
if exist(savepath,'file') && ~force_run
    return
end

basename = basenameFromBasepath(basepath);

% load needed data
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']));
load(fullfile(basepath,[basename '.ripples.events.mat']));

% make sure there are at least 5 of a each cell category available
good_to_run = check_unit_counts_per_group(cell_metrics,...
    restrict_to_brainregion,restrict_to_celltype,...
    binary_class_variable,grouping_names);
if ~good_to_run
    return
end

% restrict ripples by duration
ripples = restrict_ripples_by_duration(ripples,ripple_duration_restrict);

% get ripple spikes
ripSpk = load_spikes_and_get_ripSpk(basepath,ripples,...
    restrict_to_brainregion,restrict_to_celltype);

% restict by number of unique units
ripSpk = restrict_ripples_unique_units(ripSpk,cell_metrics,grouping_names,...
    binary_class_variable,unique_unit_num);

% get group isi distributions of each group and across groups
[A,B,AB] = calc_isi(ripSpk,cell_metrics,binary_class_variable,...
    grouping_names);

% get shuffled distributions
if shuffle
    [A_shuff,B_shuff,AB_shuff] = shuffle_isi(ripSpk,...
        cell_metrics,binary_class_variable,grouping_names);
else
    A_shuff = NaN;
    B_shuff = NaN;
    AB_shuff = NaN;
end
% save results as csv to savepath
save_results(A,B,AB,A_shuff,B_shuff,AB_shuff,grouping_names,savepath)
end

function good_to_run = check_unit_counts_per_group(cell_metrics,...
    restrict_to_brainregion,restrict_to_celltype,...
    binary_class_variable,grouping_names)

good_to_run = false;
if ~isempty(restrict_to_brainregion) && ~isempty(restrict_to_celltype)
    idx = contains(cell_metrics.brainRegion,restrict_to_brainregion) &...
        contains(cell_metrics.putativeCellType,restrict_to_celltype);
    
    if sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{1})) >= 5 &&...
            sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{2})) >= 5
        good_to_run = true;
    end
elseif ~isempty(restrict_to_brainregion) % just restrict brain region
    idx = contains(cell_metrics.brainRegion,restrict_to_brainregion);
    
    if sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{1})) >= 5 &&...
            sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{2})) >= 5
        good_to_run = true;
    end
elseif ~isempty(restrict_to_celltype) % just restrict cell type
    idx = contains(cell_metrics.putativeCellType,restrict_to_celltype);
    
    if sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{1})) >= 5 &&...
            sum(strcmp(cell_metrics.(binary_class_variable)(idx),grouping_names{2})) >= 5
        good_to_run = true;
    end
else % no restriction
    if sum(strcmp(cell_metrics.(binary_class_variable),grouping_names{1})) >= 5 &&...
            sum(strcmp(cell_metrics.(binary_class_variable),grouping_names{2})) >= 5
        good_to_run = true;
    end
end
end

function ripples = restrict_ripples_by_duration(ripples,ripple_duration_restrict)

keep = ripples.duration > ripple_duration_restrict(1) &...
    ripples.duration < ripple_duration_restrict(2);

ripples.timestamps = ripples.timestamps(keep,:);
ripples.peaks = ripples.peaks(keep,:);
ripples.amplitude = ripples.amplitude(keep,:);
ripples.frequency = ripples.frequency(keep,:);
ripples.duration = ripples.duration(keep,:);
end

function ripSpk = load_spikes_and_get_ripSpk(basepath,ripples,...
    restrict_to_brainregion,restrict_to_celltype)

% if restricting brain and cell type
if ~isempty(restrict_to_brainregion) && ~isempty(restrict_to_celltype)
    spk = importSpikes('basepath',basepath,...
        'brainRegion',restrict_to_brainregion,...
        'cellType',restrict_to_celltype);
elseif ~isempty(restrict_to_brainregion) % just restrict brain region
    spk = importSpikes('basepath',basepath,...
        'brainRegion',restrict_to_brainregion);
elseif ~isempty(restrict_to_celltype) % just restrict cell type
    spk = importSpikes('basepath',basepath,...
        'cellType',restrict_to_celltype);
else % no restriction
    spk = importSpikes('basepath',basepath);
end

% make ripSpk struct with spike times per ripple
ripSpk = getRipSpikes('basepath',basepath,'events',ripples,'spikes',spk,...
    'saveMat',false);
end

function ripSpk = restrict_ripples_unique_units(ripSpk,cell_metrics,...
    grouping_names,binary_class_variable,unique_unit_num)

for i = 1:length(ripSpk.EventRel)
    [~,Locb] = ismember(ripSpk.EventRel{i}(2,:),cell_metrics.UID);
    keep(i) = sum(strcmp(cell_metrics.(binary_class_variable)(Locb), grouping_names{1})) > unique_unit_num &&...
        sum(strcmp(cell_metrics.(binary_class_variable)(Locb), grouping_names{2})) > unique_unit_num;
end

ripSpk.EventDuration = ripSpk.EventDuration(keep);
ripSpk.UnitEventAbs = ripSpk.UnitEventAbs(:,keep);
ripSpk.UnitEventRel = ripSpk.UnitEventRel(:,keep);
ripSpk.EventAbs = ripSpk.EventAbs(keep);
ripSpk.EventRel = ripSpk.EventRel(keep);
end

function [A,B,AB] = calc_isi(ripSpk,cell_metrics,binary_class_variable,...
    grouping_names)
% iter through each rip
A=[];   B=[];   AB=[];
for e = 1:numel(ripSpk.EventRel)
    for i = 1:size(ripSpk.EventRel{e},2)-1
        for j = 1:size(ripSpk.EventRel{e},2)-1
            [A,B,AB] = get_isi(ripSpk,...
                cell_metrics,...
                j,i,e,...
                A,B,AB,...
                binary_class_variable,...
                grouping_names);
        end
    end
end
end

function [A,B,AB] = get_isi(ripSpk,cell_metrics,j,i,e,A,B,AB,...
    binary_class_variable,grouping_names)
% get_isi: main function to calc isi on sub populations

if j ~= i && ripSpk.EventRel{e}(2,i) ~= ripSpk.EventRel{e}(2,j)
    if strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)})
        
        A=cat(1,A,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
        
    elseif strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)})
        
        B=cat(1,B,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
        
    elseif (strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)})) || ...
            (strcmp(grouping_names{2},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{1},cell_metrics.(binary_class_variable){cell_metrics.UID==ripSpk.EventRel{e}(2,j)}))
        
        AB=cat(1,AB,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
    end
end
end

function [A_shuff_all,B_shuff_all,AB_shuff_all] = shuffle_isi(ripSpk,...
    cell_metrics,binary_class_variable,grouping_names)

shuff_ripSpk = ripSpk;
% -0.01 to 0.01 seconds (10ms is ~ the longest ripple cycle)
min_s = -1/100;
max_s = 1/100;

A_shuff_all = [];
B_shuff_all = [];
AB_shuff_all = [];

% 10 shuffles is enough to de-couple cycle nesting of units
for shuff_i = 1:10
    for i = 1:length(ripSpk.EventRel)
        offset = (max_s-min_s).*rand(1,size(shuff_ripSpk.EventRel{i},2)) + min_s;
        shuff_ripSpk.EventRel{i}(1,:) = ripSpk.EventRel{i}(1,:) + offset;
    end
    [A_shuff,B_shuff,AB_shuff] = calc_isi(shuff_ripSpk,...
        cell_metrics,...
        binary_class_variable,...
        grouping_names);
    A_shuff_all = [A_shuff_all;A_shuff];
    B_shuff_all = [B_shuff_all;B_shuff];
    AB_shuff_all = [AB_shuff_all;AB_shuff];
end
end

function save_results(A,B,AB,A_shuff,B_shuff,AB_shuff,grouping_names,savepath)

% nan matrix
isi = nan(max([length(A),length(B),length(AB),...
    length(A_shuff),length(B_shuff),length(AB_shuff)]),6);

% fill each column of matrix
isi(1:numel(A),1) = A;
isi(1:numel(B),2) = B;
isi(1:numel(AB),3) = AB;
isi(1:numel(A_shuff),4) = A_shuff;
isi(1:numel(B_shuff),5) = B_shuff;
isi(1:numel(AB_shuff),6) = AB_shuff;

% column names
varnames = {grouping_names{1},...
    grouping_names{2},...
    [grouping_names{1},grouping_names{2}],...
    [grouping_names{1},'_shuff'],...
    [grouping_names{2},'_shuff'],...
    [grouping_names{1},grouping_names{2},'_shuff']};

% make table
df = table(isi(:,1),isi(:,2),isi(:,3),isi(:,4),isi(:,5),isi(:,6),...
    'VariableNames',...
    varnames);

writetable(df,savepath)
end