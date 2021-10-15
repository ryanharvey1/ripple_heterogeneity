function swr_seq_pop_isi(varargin)

p = inputParser;
addParameter(p,'basepath',[]) % single basepath to run 1 session
addParameter(p,'df',[]) % data frame with df.basepath
addParameter(p,'binary_class_variable','deepSuperficial',@ischar)
addParameter(p,'grouping_names',{'Deep','Superficial'},@iscell)
addParameter(p,'restrict_to_brainregion','CA1',@ischar)
addParameter(p,'restrict_to_celltype','Pyramidal Cell',@ischar)
addParameter(p,'savepath',[])

parse(p,varargin{:})
basepaths = p.Results.basepath;
df = p.Results.df;
binary_class_variable = p.Results.binary_class_variable;
grouping_names = p.Results.grouping_names;
restrict_to_brainregion = p.Results.restrict_to_brainregion;
restrict_to_celltype = p.Results.restrict_to_celltype;
savepath = p.Results.savepath;


if isempty(df)
    df = readtable('D:\projects\ripple_heterogeneity\sessions.csv');
    basepaths = unique(df.basepath);
end

WaitMessage = parfor_wait(length(basepaths));
parfor i = 1:length(basepaths)
    basepath = basepaths{i};
    disp(basepath)
    
    file_str = strsplit(basepath,filesep);
    savepath_ = fullfile(savepath,[file_str{end-1},'_',file_str{end},'.csv']);
    if exist(savepath_,'file')
        continue
    end
    
    main(basepath,...
        binary_class_variable,...
        grouping_names,...
        restrict_to_brainregion,...
        restrict_to_celltype,...
        savepath_)
    
    WaitMessage.Send;
end
WaitMessage.Destroy;
end

function main(basepath,binary_class_variable,grouping_names,...
    restrict_to_brainregion,restrict_to_celltype,savepath)

basename = basenameFromBasepath(basepath);

% load needed data
load(fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']));
load(fullfile(basepath,[basename '.ripples.events.mat']));

ripSpk = load_spikes_and_get_ripSpk(basepath,ripples,...
    restrict_to_brainregion,restrict_to_celltype);

% iter through each rip
A=[];   B=[];   AB=[];
for e = 1:numel(ripSpk.EventRel)
    for i = 1:length(ripSpk.EventRel{e})-1
        for j = 1:length(ripSpk.EventRel{e})-1
            [A,B,AB] = get_isi(ripSpk,...
                cell_metrics,...
                j,i,e,...
                A,B,AB,...
                binary_class_variable,...
                grouping_names);
        end
    end
end
save_results(A,B,AB,grouping_names,savepath)
end

function save_results(A,B,AB,grouping_names,savepath)

meanISI = nan(max([length(A),length(B),length(AB)]),3);
meanISI(1:numel(A),1) = A;
meanISI(1:numel(B),2) = B;
meanISI(1:numel(AB),3) = AB;

df = table(meanISI(:,1),meanISI(:,2),meanISI(:,3),...
    'VariableNames',...
    {grouping_names{1},grouping_names{2},[grouping_names{1},grouping_names{2}]});

writetable(df,savepath)
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
ripSpk = getRipSpikes('basepath',basepath,'events',ripples,'spikes',spk,'saveMat',false);

end

function [A,B,AB] = get_isi(ripSpk,cell_metrics,j,i,e,A,B,AB,...
    binary_class_variable,grouping_names)
% get_isi: main function to calc isi on sub populations

if j ~= i && ripSpk.EventRel{e}(2,i) ~= ripSpk.EventRel{e}(2,j)
    if strcmp(grouping_names{1},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{1},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,j)})
        
        A=cat(1,A,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
        
    elseif strcmp(grouping_names{2},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{2},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,j)})
        
        B=cat(1,B,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
        
    elseif (strcmp(grouping_names{1},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{2},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,j)})) || ...
            (strcmp(grouping_names{2},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,i)}) && ...
            strcmp(grouping_names{1},cell_metrics.(binary_class_variable){ripSpk.EventRel{e}(2,j)}))
        
        AB=cat(1,AB,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
    end
end
end