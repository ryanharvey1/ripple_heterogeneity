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

%% count cells 
nCA1pyr = 0; ses = 0; nCells = 0;
for a = 1:numel(animal)
    for d = 1:numel(day{a})
        if strncmp('OML',animal{a},3)
            cd([dataDir2 animal{a} '\' day{a}{d}]); 
        elseif strncmp('Wmaze',animal{a},5)
            cd([dataDir3 animal{a} '\' day{a}{d}]); 
        else
            cd([dataDir1 animal{a} '\' day{a}{d}]);
        end
        ses = ses+1;
        basename = bz_BasenameFromBasepath(pwd);
        if exist([basename '.cell_metrics.cellinfo.mat'],'file')
            load([basename '.cell_metrics.cellinfo.mat']);
        else
            disp([basename '.cell_metrics.cellinfo.mat',' does not exist'])
            continue
        end
       
        ntemp=0;
        for i = 1:numel(cell_metrics.putativeCellType) 
            if strcmp(cell_metrics.putativeCellType{i},'Pyramidal Cell') && ...
               strcmp(cell_metrics.brainRegion{i},'CA1')
               ntemp= ntemp+1;
            end
        end
        nCA1pyrSes(ses,1) = ntemp;
        nCA1pyr = nCA1pyr + ntemp;
        nCells = nCells + numel(cell_metrics.putativeCellType);
        clear ntemp cell_metrics;
    end
end

%% pool data v1
sesPre=0; sesTask=0; sesPost=0;
for a = 1:numel(animal)%13
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
            disp([basename '.SWRunitMetrics.mat',' does not exist'])
            continue
        end
        if isfield(SWRunitMetrics,'pre') && ~isempty(SWRunitMetrics.pre)
           sesPre = sesPre+1;
           SWRunitMetrics.pre.cellType = cell_metrics.putativeCellType;
           SWRunitMetrics.pre.region = cell_metrics.brainRegion;
           SWRunitMetrics.pre.CA1depth = cell_metrics.CA1depth;
           % other cell metrics of interest: FR, burst, etc. 
           unitPre{sesPre} = SWRunitMetrics.pre;
        end
        if isfield(SWRunitMetrics,'task') && ~isempty(SWRunitMetrics.task)
           sesTask = sesTask+1;
           SWRunitMetrics.task.cellType = cell_metrics.putativeCellType;
           SWRunitMetrics.task.region = cell_metrics.brainRegion;
           SWRunitMetrics.task.CA1depth = cell_metrics.CA1depth;           
           unitTask{sesTask} = SWRunitMetrics.task;
        end        
        if isfield(SWRunitMetrics,'post') && ~isempty(SWRunitMetrics.post)
           sesPost = sesPost+1;
           SWRunitMetrics.post.cellType = cell_metrics.putativeCellType;
           SWRunitMetrics.post.region = cell_metrics.brainRegion;
           SWRunitMetrics.post.CA1depth = cell_metrics.CA1depth;           
           unitPost{sesPost} = SWRunitMetrics.post;
        end        
        
        clear SWRunitMetrics
    end
end

unitPre = cell2mat(unitPre);
unitTask = cell2mat(unitTask);
unitPost = cell2mat(unitPost);

allUnitPre = bz_CollapseStruct(unitPre,'match');
allUnitTask = bz_CollapseStruct(unitTask,'match');
allUnitPost = bz_CollapseStruct(unitPost,'match');

%% clasify cell types v1
% 1- Deep, middle, sup
for i = 1:numel(allUnitPre.particip)
    if strcmp(allUnitPre.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitPre.region{i},'CA1') && ...   
       allUnitPre.CA1depth(i) > 0 
       allUnitPre.ID(i,1) = 1;
    elseif strcmp(allUnitPre.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitPre.region{i},'CA1') && ...   
       allUnitPre.CA1depth(i) == 0 
       allUnitPre.ID(i,1) = 2;
    elseif strcmp(allUnitPre.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitPre.region{i},'CA1') && ...   
       allUnitPre.CA1depth(i) < 0 
       allUnitPre.ID(i,1) = 3;
    else
        allUnitPre.ID(i,1) = 0;
    end
end

for i = 1:numel(allUnitTask.particip)
    if strcmp(allUnitTask.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitTask.region{i},'CA1') && ...   
       allUnitTask.CA1depth(i) > 0 
       allUnitTask.ID(i,1) = 1;
    elseif strcmp(allUnitTask.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitTask.region{i},'CA1') && ...   
       allUnitTask.CA1depth(i) == 0 
       allUnitTask.ID(i,1) = 2;
    elseif strcmp(allUnitTask.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitTask.region{i},'CA1') && ...   
       allUnitTask.CA1depth(i) < 0 
       allUnitTask.ID(i,1) = 3;
    else
        allUnitTask.ID(i,1) = 0;
    end
end

for i = 1:numel(allUnitPost.particip)
    if strcmp(allUnitPost.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitPost.region{i},'CA1') && ...   
       allUnitPost.CA1depth(i) > 0 
       allUnitPost.ID(i,1) = 1;
    elseif strcmp(allUnitPost.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitPost.region{i},'CA1') && ...   
       allUnitPost.CA1depth(i) == 0 
       allUnitPost.ID(i,1) = 2;
    elseif strcmp(allUnitPost.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitPost.region{i},'CA1') && ...   
       allUnitPost.CA1depth(i) < 0 
       allUnitPost.ID(i,1) = 3;
    else
        allUnitPost.ID(i,1) = 0;
    end
end

% 2 - REM shifting / non shifting 

%% pool data v2 
sesAll = 0;
for a = 1:numel(animal)
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
            disp([basename '.SWRunitMetrics.mat',' does not exist'])
            continue
        end
        sesAll = sesAll+1;
           unitAll{sesAll}.cellType = cell_metrics.putativeCellType;
           unitAll{sesAll}.region = cell_metrics.brainRegion;
           unitAll{sesAll}.CA1depth = cell_metrics.CA1depth;        
        
        if isfield(SWRunitMetrics,'pre') && ~isempty(SWRunitMetrics.pre)
           unitAll{sesAll}.particip(:,1) = SWRunitMetrics.pre.particip;
           unitAll{sesAll}.FRall(:,1) = SWRunitMetrics.pre.particip;
        else
           unitAll{sesAll}.particip(1:numel(cell_metrics.brainRegion),1) = NaN;
           unitAll{sesAll}.FRall(1:numel(cell_metrics.brainRegion),1) = NaN;          
        end
        if isfield(SWRunitMetrics,'task') && ~isempty(SWRunitMetrics.task)
           unitAll{sesAll}.particip(:,2) = SWRunitMetrics.task.particip;
           unitAll{sesAll}.FRall(:,2) = SWRunitMetrics.task.particip;
        else
           unitAll{sesAll}.particip(1:numel(cell_metrics.brainRegion),2) = NaN;
           unitAll{sesAll}.FRall(1:numel(cell_metrics.brainRegion),2) = NaN;          
        end
        if isfield(SWRunitMetrics,'post') && ~isempty(SWRunitMetrics.post)
           unitAll{sesAll}.particip(:,3) = SWRunitMetrics.post.particip;
           unitAll{sesAll}.FRall(:,3) = SWRunitMetrics.post.particip;
        else
           unitAll{sesAll}.particip(1:numel(cell_metrics.brainRegion),3) = NaN;
           unitAll{sesAll}.FRall(1:numel(cell_metrics.brainRegion),3) = NaN;          
        end
        
    end
end

unitAll = cell2mat(unitAll);
allUnitAll = bz_CollapseStruct(unitAll,'match');

% clasify cell types v2
for i = 1:length(allUnitAll.particip)
    if strcmp(allUnitAll.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitAll.region{i},'CA1') && ...   
       allUnitAll.CA1depth(i) > 0 
       allUnitAll.ID(i,1) = 1;
    elseif strcmp(allUnitAll.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitAll.region{i},'CA1') && ...   
       allUnitAll.CA1depth(i) == 0 
       allUnitAll.ID(i,1) = 2;
    elseif strcmp(allUnitAll.cellType{i},'Pyramidal Cell') && ...
       strcmp(allUnitAll.region{i},'CA1') && ...   
       allUnitAll.CA1depth(i) < 0 
       allUnitAll.ID(i,1) = 3;
    else
        allUnitAll.ID(i,1) = 0;
    end
end

%% plot
binrange = 0:0.05:1; f = 3;
figure;
plot(binrange,smooth(histc(allUnitPre.particip(allUnitPre.ID==1),binrange)/sum(histc(allUnitPre.particip(allUnitPre.ID==1),binrange)),f,'sgolay'),'r');hold on;
plot(binrange,smooth(histc(allUnitPre.particip(allUnitPre.ID==2),binrange)/sum(histc(allUnitPre.particip(allUnitPre.ID==2),binrange)),f,'sgolay'),'k');hold on;
plot(binrange,smooth(histc(allUnitPre.particip(allUnitPre.ID==3),binrange)/sum(histc(allUnitPre.particip(allUnitPre.ID==3),binrange)),f,'sgolay'),'b');hold on;
plot([median(allUnitPre.particip(allUnitPre.ID==1)) median(allUnitPre.particip(allUnitPre.ID==1))],ylim,'r'); hold on;
plot([median(allUnitPre.particip(allUnitPre.ID==2)) median(allUnitPre.particip(allUnitPre.ID==2))],ylim,'k'); hold on;
plot([median(allUnitPre.particip(allUnitPre.ID==3)) median(allUnitPre.particip(allUnitPre.ID==3))],ylim,'b'); hold on;
ylim([0 Inf]);xlim([0 1]); xlabel('particip. prob'); ylabel('frac. of cells');legend({'deep','mid','sup'});
p=ranksum(allUnitPre.particip(allUnitPre.ID==1),allUnitPre.particip(allUnitPre.ID==3));title(['p = ' num2str(p,2)]);
%set(gca,'XScale','log')

figure;
plot(binrange,smooth(histc(allUnitPost.particip(allUnitPost.ID==1),binrange)/sum(histc(allUnitPost.particip(allUnitPost.ID==1),binrange)),f,'sgolay'),'r');hold on;
plot(binrange,smooth(histc(allUnitPost.particip(allUnitPost.ID==2),binrange)/sum(histc(allUnitPost.particip(allUnitPost.ID==2),binrange)),f,'sgolay'),'k');hold on;
plot(binrange,smooth(histc(allUnitPost.particip(allUnitPost.ID==3),binrange)/sum(histc(allUnitPost.particip(allUnitPost.ID==3),binrange)),f,'sgolay'),'b');hold on;
plot([median(allUnitPost.particip(allUnitPost.ID==1)) median(allUnitPost.particip(allUnitPost.ID==1))],ylim,'r'); hold on;
plot([median(allUnitPost.particip(allUnitPost.ID==2)) median(allUnitPost.particip(allUnitPost.ID==2))],ylim,'k'); hold on;
plot([median(allUnitPost.particip(allUnitPost.ID==3)) median(allUnitPost.particip(allUnitPost.ID==3))],ylim,'b'); hold on;
ylim([0 Inf]);xlim([0 1]); xlabel('particip. prob'); ylabel('frac. of cells');legend({'deep','mid','sup'});
p=ranksum(allUnitPost.particip(allUnitPost.ID==1),allUnitPost.particip(allUnitPost.ID==3));title(['p = ' num2str(p,2)]);
%set(gca,'XScale','log')

%%
binrange = 0:0.05:1; f = 3;
figure;
subplot(1,2,1);
plot(binrange,smooth(histc(allUnitAll.particip(allUnitAll.ID>0,1),binrange)/sum(histc(allUnitAll.particip(allUnitAll.ID>0,1),binrange)),f,'sgolay'),'k');hold on;
ylim([0 Inf]);xlim([0 1]); xlabel('particip. prob'); ylabel('frac. of cells');legend({'all CA1pyr'});
subplot(1,2,2);
scatter(allUnitAll.particip(allUnitAll.ID>0,1),allUnitAll.particip(allUnitAll.ID>0,2),'.');hold on;
[r p]=corr(allUnitAll.particip(allUnitAll.ID>0,1),allUnitAll.particip(allUnitAll.ID>0,2),'rows','complete');
title(['r = ' num2str(r,2) ' / p = ' num2str(p,2)]);
plot([0 1],[0 1],'--r'); xlabel('particp pre');ylabel('particp task');
subplot(1,3,3);
scatter(allUnitAll.particip(allUnitAll.ID>0,1),allUnitAll.particip(allUnitAll.ID>0,3),'.');hold on;
[r p]=corr(allUnitAll.particip(allUnitAll.ID>0,1),allUnitAll.particip(allUnitAll.ID>0,3),'rows','complete');
title(['r = ' num2str(r,2) ' / p = ' num2str(p,2)]);
plot([0 1],[0 1],'--r'); xlabel('particp pre');ylabel('particp post');


%% write to table

varNames = {'particip','FRall','FRparticip','nSpkAll','nSpkParticip',...
    'cellType','region','CA1depth','ID'};

df_pre = table(allUnitPre.particip,...
    allUnitPre.FRall,...
    allUnitPre.FRparticip,...
    allUnitPre.nSpkAll,...
    allUnitPre.nSpkParticip,...
    allUnitPre.cellType',...
    allUnitPre.region',...
    allUnitPre.CA1depth,...
    allUnitPre.ID,...
    'VariableNames',varNames);
df_pre.epoch(:) = {'pre'};

df_task = table(allUnitTask.particip,...
    allUnitTask.FRall,...
    allUnitTask.FRparticip,...
    allUnitTask.nSpkAll,...
    allUnitTask.nSpkParticip,...
    allUnitTask.cellType',...
    allUnitTask.region',...
    allUnitTask.CA1depth,...
    allUnitTask.ID,...
    'VariableNames',varNames);
df_task.epoch(:) = {'task'};

df_post = table(allUnitPost.particip,...
    allUnitPost.FRall,...
    allUnitPost.FRparticip,...
    allUnitPost.nSpkAll,...
    allUnitPost.nSpkParticip,...
    allUnitPost.cellType',...
    allUnitPost.region',...
    allUnitPost.CA1depth,...
    allUnitPost.ID,...
    'VariableNames',varNames);
df_post.epoch(:) = {'post'};

df = [df_pre;df_task;df_post];

% assign ID to real name
df.ca1_layer = df.cellType; % placeholder
df.ca1_layer(df.ID == 1) = {'deep'};
df.ca1_layer(df.ID == 2) = {'mid'};
df.ca1_layer(df.ID == 1) = {'sup'};
df.ca1_layer(df.ID == 0) = {'unknown'};

writetable(df,'D:\projects\ripple_heterogeneity\df.csv')
% varNames = {'cellType','region','CA1depth','particip','FRall','ID'};
%  
% df_all_unit = table(allUnitAll.cellType',...
%     allUnitAll.region',...
%     allUnitAll.CA1depth,...
%     allUnitAll.particip(:,1),...
%     allUnitAll.FRall(:,1),...
%     allUnitAll.ID,...
%     'VariableNames',varNames);
% 
% add different tasks (pre task post)
% df_all_unit.epoch(:) = {'pre'};
% 
% df_all_unit_task = table(allUnitAll.cellType',...
%     allUnitAll.region',...
%     allUnitAll.CA1depth,...
%     allUnitAll.particip(:,2),...
%     allUnitAll.FRall(:,2),...
%     allUnitAll.ID,...
%     'VariableNames',varNames);
% 
% df_all_unit_post = table(allUnitAll.cellType',...
%     allUnitAll.region',...
%     allUnitAll.CA1depth,...
%     allUnitAll.particip(:,3),...
%     allUnitAll.FRall(:,3),...
%     allUnitAll.ID,...
%     'VariableNames',varNames);


