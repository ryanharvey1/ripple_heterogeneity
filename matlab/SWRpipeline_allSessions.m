%% all control sessions processed for SWRpipeline analysis
%  FARNAZ: add here Grossmark sessions and OML18/19
clearvars;
dataDir1 = 'A:\Data\';
dataDir2 = 'A:\OptoMECLEC\';
dataDir3 = 'A:\ORproject\';
animal = {'AB1','AB3','AB4','AYA4','AYA6','AYA7','AYA9','AYA10','OML5','OML3','OML7','OML8','OML10',...
    'OML18','OML19','Wmaze2\OR15','Wmaze2\OR18','Wmaze3\OR22','Wmaze3\OR21','Wmaze3\OR23',...
    'GrosmarkAD\Cicero','GrosmarkAD\Buddy','GrosmarkAD\Achilles','GrosmarkAD\Gatsby'};
day = {{'day1'},{'AB3_38_41','AB3_42_46','AB3_47_49','AB3_50_51','AB3_55_57','AB3_58_59','AB3_60'},...
    {'day03','day07','day08','day09','day11'},{'day150726','day150728','day150804'},...
    {'day17','day19','day20'},{'day19','day20','day22','day24','day25','day27','day30'},...
    {'day12','day15','day16','day17','day20'},{'day25','day27','day31','day32','day34'},...
    {'day4','day5','day6','day8','day9','day12','day13','day20','day21'},...
    {'day11','day12','day13','day14','day16','day17'},{'day6','day7','day9'},...
    {'day4','day5','day6','day7','day8','day9','day17','day20'},{'day5','day7','day9'},...
    {'day1','day2','day4','day5'},{'day2','day3'},...
    {'day1','day2','day3','day4','day10'},{'day1','day2','day3'},{'day1','day4','day3','day5'},...
    {'day2','day4'},{'day1','day5'},...
    {'Cicero_09012014','Cicero_09102014','Cicero_09172014'},{'Buddy_06272013'},{'Achilles_10252013','Achilles_11012013'},{'Gatsby_08282013','Gatsby_08022013'}};
opto =    {NaN,[NaN NaN NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],[NaN NaN NaN],...
    [NaN NaN NaN],[NaN NaN NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],...
    [1 0 1 1 0 8 3 0 1],[1 0 1 1 3 3],[1 0 1],[0 2 0 2 0 2 2 2],[2 2 2],[1 1 1 1],[1 1]...
    [4 4 4 4 4],[4 4 4],[4 4 4 4],[4 4],[4 4],[NaN NaN NaN],NaN,[NaN NaN],[NaN NaN]};
% 0=sham, 1=MEC silencing; 2=LEC silencing; 3=prb stm; 4=SWR prolong; 8=problem
task =    {1,[1 1 4 4 4 1 4],[2 2 2 2 2],[1 1 1],[1 1 1],[1 1 2 2 2 1 2],[1 2 2 1 2],[2 2 2 2 2],...
    [2 2 2 2 2 2 2 2 2],[2 2 2 2 2 2],[2 2 2],[2 2 2 2 2 2 2 2],[2 2 2],[1 1 1 1],[1 1]...
    [3 3 3 3 3],[3 3 3],[3 3 3 3],[3 3],[3 3],[1 1 1],1,[1 1],[1 1]};
% 1=linear, 2=cheesboard, 3=Wmaze, 4=Tmaze
novel =   {0,[0 0 0 0 0 0 0],[0 0 0 0 0],[0 0 0],[0 0 0],[0 0 0 0 0 0 0],[0 0 0 0 0],[0 0 0 0 0],...
    [0 0 0 0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0],[0 0 0 0 0 0 0 0],[0 0 0],[0 0 0 0],[0 0]...
    [0 0 0 0 0],[0 0 0],[0 0 0 0],[0 0],[0 0],[1 1 1],1,[1 1],[1 1]};
% 0=familiar env, 1=novel env.   THIS NEEDS TO BE CHECKED, THERE ARE SOME NOVEL SES
regions = {1,[1 1 1 1 1 1 1],[1 1 1 1 1],[3 3 3],[1 1 1],[2 2 2 2 2 2 2],[2 2 2 2 2],[3 3 3 3 3],...
    [1 1 1 1 1 1 1 1 1],[1 1 1 1 1 1],[1 1 1],[1 1 1 1 1 1 1 1],[1 1 1],[1 1 1 1],[1 1]...
    [1 1 1 1 1],[4 4 4],[1 1 1 1],[1 1],[1 1],[1 1 1],1,[1 1],[1 1]};
% 1=only HP, 2=also MEC, 3=also LEC, 4=PFC

% PROBLEMS: OML8day17 (States)
% PENDING: Gabi, Kenji, Shige
% pyr/int clasific OML and OR


%% pool data v1
sesPre=0; sesTask=0; sesPost=0;
need_attention = [];
check_cell_metrics = [];
for a = 1:numel(animal)%13
    disp(animal{a})
    for d = 1:numel(day{a})
        if strncmp('OML',animal{a},3)
            cd([dataDir2 animal{a} '\' day{a}{d}]);
        elseif strncmp('Wmaze',animal{a},5)
            cd([dataDir3 animal{a} '\' day{a}{d}]);
        else
            cd([dataDir1 animal{a} '\' day{a}{d}]);
        end
        basename = bz_BasenameFromBasepath(pwd);
        if exist([basename '.SWRunitMetrics.mat'],'file')
            load([basename '.SWRunitMetrics.mat']);
            load([basename '.cell_metrics.cellinfo.mat']);
        else
            need_attention = [need_attention;{[animal{a},'_',day{a}{d}]}];
            continue
        end
        
        if isfield(SWRunitMetrics,'pre') && ~isempty(SWRunitMetrics.pre)
            sesPre = sesPre+1;
            [SWRunitMetrics,needs_cell_metrics] = extract_from_metrics(cell_metrics,...
                SWRunitMetrics,'pre',a,d,animal,day,opto,task,novel,regions);
            
            % other cell metrics of interest: FR, burst, etc.
            unitPre{sesPre} = SWRunitMetrics.pre;
        end
        if isfield(SWRunitMetrics,'task') && ~isempty(SWRunitMetrics.task)
            sesTask = sesTask+1;
            
            [SWRunitMetrics,needs_cell_metrics] = extract_from_metrics(cell_metrics,...
                SWRunitMetrics,'task',a,d,animal,day,opto,task,novel,regions);
            
            unitTask{sesTask} = SWRunitMetrics.task;
        end
        if isfield(SWRunitMetrics,'post') && ~isempty(SWRunitMetrics.post)
            sesPost = sesPost+1;
            
            [SWRunitMetrics,needs_cell_metrics] = extract_from_metrics(cell_metrics,...
                SWRunitMetrics,'post',a,d,animal,day,opto,task,novel,regions);
            
            unitPost{sesPost} = SWRunitMetrics.post;
        end
        
        check_cell_metrics = [check_cell_metrics;{needs_cell_metrics}];
        
        clear SWRunitMetrics
    end
end

disp('needs cell metrics update')
disp(check_cell_metrics(~cellfun('isempty',check_cell_metrics)))

disp('needs SWRunitMetrics update')
disp(need_attention)

unitPre = cell2mat(unitPre);
unitTask = cell2mat(unitTask);
unitPost = cell2mat(unitPost);

allUnitPre = bz_CollapseStruct(unitPre,'match');
allUnitTask = bz_CollapseStruct(unitTask,'match');
allUnitPost = bz_CollapseStruct(unitPost,'match');

%% clasify cell types v1
% 1- Deep, middle, sup
allUnitPre.ID(allUnitPre.CA1depth>0 & strcmp(allUnitPre.region,'CA1')',1) = 1;
allUnitPre.ID(allUnitPre.CA1depth==0 & strcmp(allUnitPre.region,'CA1')',1) = 2;
allUnitPre.ID(allUnitPre.CA1depth<0 & strcmp(allUnitPre.region,'CA1')',1) = 3;

allUnitTask.ID(allUnitTask.CA1depth>0 & strcmp(allUnitTask.region,'CA1')',1) = 1;
allUnitTask.ID(allUnitTask.CA1depth==0 & strcmp(allUnitTask.region,'CA1')',1) = 2;
allUnitTask.ID(allUnitTask.CA1depth<0 & strcmp(allUnitTask.region,'CA1')',1) = 3;

allUnitPost.ID(allUnitPost.CA1depth>0 & strcmp(allUnitPost.region,'CA1')',1) = 1;
allUnitPost.ID(allUnitPost.CA1depth==0 & strcmp(allUnitPost.region,'CA1')',1) = 2;
allUnitPost.ID(allUnitPost.CA1depth<0 & strcmp(allUnitPost.region,'CA1')',1) = 3;


%% write to table

varNames = {'particip','FRall','FRparticip','nSpkAll','nSpkParticip',...
    'cellType','region','CA1depth','UID','burstIndex_Royer2012',...
    'animal','day','opto','task','novel','ID'};

df_pre = make_table(allUnitPre,varNames);
df_pre.epoch(:) = {'pre'};

df_task = make_table(allUnitTask,varNames);
df_task.epoch(:) = {'task'};

df_post = make_table(allUnitPost,varNames);
df_post.epoch(:) = {'post'};

df = [df_pre;df_task;df_post];


%% assign labels to dummy coded values

% assign ID to real name
df.ca1_layer = df.cellType; % placeholder
df.ca1_layer(df.ID == 0) = {'unknown'};
df.ca1_layer(df.ID == 1) = {'deep'};
df.ca1_layer(df.ID == 2) = {'mid'};
df.ca1_layer(df.ID == 3) = {'sup'};
df.ID = [];

% 0=sham, 1=MEC silencing; 2=LEC silencing; 3=prb stm; 4=SWR prolong; 8=problem
temp_var = repmat({'none'},length(df.opto),1);
temp_var(df.opto == 0) = {'sham'};
temp_var(df.opto == 1) = {'mec_silencing'};
temp_var(df.opto == 2) = {'lec_silencing'};
temp_var(df.opto == 3) = {'prb_stm'};
temp_var(df.opto == 4) = {'swr_prolong'};
temp_var(df.opto == 8) = {'problem'};
df.opto = temp_var;

% 1=linear, 2=cheesboard, 3=Wmaze, 4=Tmaze
temp_var = repmat({'unknown'},length(df.task),1);
temp_var(df.task == 1) = {'linear'};
temp_var(df.task == 2) = {'cheesboard'};
temp_var(df.task == 3) = {'w_maze'};
temp_var(df.task == 4) = {'t_maze'};
df.task = temp_var;

% 0=familiar env, 1=novel env.
temp_var = repmat({'unknown'},length(df.novel),1);
temp_var(df.novel == 0) = {'familiar_env'};
temp_var(df.novel == 1) = {'novel_env'};
df.novel = temp_var;

%% write to csv for further processing
writetable(df,'D:\projects\ripple_heterogeneity\df.csv')





%% local functions


function [SWRunitMetrics,need_attention] = extract_from_metrics(cell_metrics,...
    SWRunitMetrics,epoch,a,d,animal,day,opto,task,novel,regions)

need_attention = [];

SWRunitMetrics.(epoch).cellType = cell_metrics.putativeCellType;
SWRunitMetrics.(epoch).region = cell_metrics.brainRegion;
SWRunitMetrics.(epoch).CA1depth = cell_metrics.CA1depth;
if isfield(cell_metrics,'UID')
    SWRunitMetrics.(epoch).UID = cell_metrics.UID(:);
    SWRunitMetrics.(epoch).burstIndex_Royer2012 = cell_metrics.burstIndex_Royer2012(:);
else
    SWRunitMetrics.(epoch).UID = nan(length(cell_metrics.CA1depth),1);
    SWRunitMetrics.(epoch).burstIndex_Royer2012 = nan(length(cell_metrics.CA1depth),1);
    need_attention = [animal{a},'_',day{a}{d}];
end
% add ID and metadata
SWRunitMetrics.(epoch).animal = repmat(animal(a),length(SWRunitMetrics.(epoch).cellType),1);
SWRunitMetrics.(epoch).day = repmat(day{a}(d),length(SWRunitMetrics.(epoch).cellType),1);
SWRunitMetrics.(epoch).opto = repmat(opto{a}(d),length(SWRunitMetrics.(epoch).cellType),1);
SWRunitMetrics.(epoch).task = repmat(task{a}(d),length(SWRunitMetrics.(epoch).cellType),1);
SWRunitMetrics.(epoch).novel = repmat(novel{a}(d),length(SWRunitMetrics.(epoch).cellType),1);
SWRunitMetrics.(epoch).regions = repmat(regions{a}(d),length(SWRunitMetrics.(epoch).cellType),1);
end

function df = make_table(structure,varNames)
df = table();
for v = varNames
    df.(v{1}) = structure.(v{1})(:);
end
end
