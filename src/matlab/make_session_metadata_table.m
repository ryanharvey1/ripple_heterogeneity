
% make_session_metadata_table
animal = {'AB1','AB3','AB4','AYA4','AYA6','AYA7','AYA9','AYA10','OML5','OML3','OML7','OML8','OML10',...
    'OML18','OML19','Wmaze2\OR15','Wmaze2\OR18','Wmaze3\OR22','Wmaze3\OR21','Wmaze3\OR23',...
    'GrosmarkAD\Cicero','GrosmarkAD\Buddy','GrosmarkAD\Achilles','GrosmarkAD\Gatsby',...
    'Kenji'};
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
    {'Cicero_09012014','Cicero_09102014','Cicero_09172014'},{'Buddy_06272013'},...
    {'Achilles_10252013','Achilles_11012013'},{'Gatsby_08282013','Gatsby_08022013'},...
    {'ec013.895_902'}};
opto =    {NaN,[NaN NaN NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],[NaN NaN NaN],...
    [NaN NaN NaN],[NaN NaN NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],...
    [1 0 1 1 0 8 3 0 1],[1 0 1 1 3 3],[1 0 1],[0 2 0 2 0 2 2 2],[2 2 2],[1 1 1 1],[1 1]...
    [4 4 4 4 4],[4 4 4],[4 4 4 4],[4 4],[4 4],[NaN NaN NaN],NaN,[NaN NaN],[NaN NaN],NaN};
% 0=sham, 1=MEC silencing; 2=LEC silencing; 3=prb stm; 4=SWR prolong; 8=problem
task =    {1,[1 1 4 4 4 1 4],[2 2 2 2 2],[1 1 1],[1 1 1],[1 1 2 2 2 1 2],[1 2 2 1 2],[2 2 2 2 2],...
    [2 2 2 2 2 2 2 2 2],[2 2 2 2 2 2],[2 2 2],[2 2 2 2 2 2 2 2],[2 2 2],[1 1 1 1],[1 1]...
    [3 3 3 3 3],[3 3 3],[3 3 3 3],[3 3],[3 3],[1 1 1],1,[1 1],[1 1],5};
% 1=linear, 2=cheesboard, 3=Wmaze, 4=Tmaze, 5=Mwheel
novel =   {0,[0 0 0 0 0 0 0],[0 0 0 0 0],[0 0 0],[0 0 0],[0 0 0 0 0 0 0],[0 0 0 0 0],[0 0 0 0 0],...
    [0 0 0 0 0 0 0 0 0],[0 0 0 0 0 0],[0 0 0],[0 0 0 0 0 0 0 0],[0 0 0],[0 0 0 0],[0 0]...
    [0 0 0 0 0],[0 0 0],[0 0 0 0],[0 0],[0 0],[1 1 1],1,[1 1],[1 1],0};
% 0=familiar env, 1=novel env.   THIS NEEDS TO BE CHECKED, THERE ARE SOME NOVEL SES
regions = {1,[1 1 1 1 1 1 1],[1 1 1 1 1],[3 3 3],[1 1 1],[2 2 2 2 2 2 2],[2 2 2 2 2],[3 3 3 3 3],...
    [1 1 1 1 1 1 1 1 1],[1 1 1 1 1 1],[1 1 1],[1 1 1 1 1 1 1 1],[1 1 1],[1 1 1 1],[1 1]...
    [1 1 1 1 1],[4 4 4],[1 1 1 1],[1 1],[1 1],[1 1 1],1,[1 1],[1 1],NaN};
% 1=only HP, 2=also MEC, 3=also LEC, 4=PFC

% unpack
i = 1;
for a = 1:length(animal)
    for d = 1:length(day{a})
        df_animal{i} = animal{a};
        df_day{i} = day{a}{d};
        df_opto(i) = opto{a}(d);
        df_task(i) = task{a}(d);
        df_novel(i) = novel{a}(d);
        df_regions(i) = regions{a}(d);
        i = i+1;
    end
end
df = table();
df.animal = df_animal';
df.day = df_day';
df.opto = df_opto';
df.task = df_task';
df.novel = df_novel';
df.region = df_regions';

% 0=sham, 1=MEC silencing; 2=LEC silencing; 3=prb stm; 4=SWR prolong; 8=problem
opto_temp = repmat({'none'},length(df.opto),1);
opto_temp(df.opto == 0) = {'sham'};
opto_temp(df.opto == 1) = {'mec_silencing'};
opto_temp(df.opto == 2) = {'lec_silencing'};
opto_temp(df.opto == 3) = {'prb_stm'};
opto_temp(df.opto == 4) = {'swr_prolong'};
opto_temp(df.opto == 8) = {'problem'};
df.opto = opto_temp;

% 1=linear, 2=cheesboard, 3=Wmaze, 4=Tmaze, 5=Mwheel
task_temp = repmat({'unknown'},length(df.task),1);
task_temp(df.task == 1) = {'linear'};
task_temp(df.task == 2) = {'cheesboard'};
task_temp(df.task == 3) = {'w_maze'};
task_temp(df.task == 4) = {'t_maze'};
task_temp(df.task == 5) = {'m_wheel'};
df.task = task_temp;

% 0=familiar env, 1=novel env.   THIS NEEDS TO BE CHECKED, THERE ARE SOME NOVEL SES
novel_temp = repmat({'unknown'},length(df.novel),1);
novel_temp(df.novel == 0) = {'familiar'};
novel_temp(df.novel == 1) = {'novel'};
df.novel = novel_temp;

% 1=only HP, 2=also MEC, 3=also LEC, 4=PFC
region_temp = repmat({'unknown'},length(df.region),1);
region_temp(df.region == 1) = {'hp'};
region_temp(df.region == 2) = {'hp_mec'};
region_temp(df.region == 3) = {'hp_lec'};
region_temp(df.region == 4) = {'pfc'};
df.region = region_temp;

df(contains(df.opto,'mec'),:)
% df(df.task==1 & df.opto == 1,:)

