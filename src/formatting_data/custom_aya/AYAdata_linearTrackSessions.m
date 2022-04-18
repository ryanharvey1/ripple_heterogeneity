%% Linear track sesssions for SWR analysis 
% all these sesssions have linear track and sleep (ideally pre and post
% sleep but not always). 

dataDir = 'A:\Data\';
animal = {'AB1','AB3','AYA4','AYA6','AYA7','AYA9'};
day = {{'day1'},{'AB3_38_41','AB3_42_46','AB3_58_59'},{'day150726','day150728','day150804'},...
       {'day17','day19','day20'},{'day19','day20','day27'},{'day12','day17'}};

%% check all files are present
for a = 1:numel(animal)
    for d = 1:numel(day{a})
        cd([dataDir animal{a} '\' day{a}{d}]);
        basename = bz_BasenameFromBasepath(pwd);
        
        if FileExists([basename '.cell_metrics.cellinfo.mat']) && ...
           FileExists([basename '.SWRunitMetrics.mat']) && FileExists('posTrials.mat') && ...
           FileExists([basename '.thetaPhaseStates.mat']) && FileExists([basename '.swrCh.mat'])
           disp([animal{a} '\' day{a}{d} '   CHECK'])
        else
           disp([animal{a} '\' day{a}{d} '   MISSING FILE'])
        end
    end
end

%% Count n CA1pyr 
nCA1pyr = 0; ses = 0;
for a = 1:numel(animal)
    for d = 1:numel(day{a})
        cd([dataDir animal{a} '\' day{a}{d}]);
        ses = ses+1;
        basename = bz_BasenameFromBasepath(pwd);
        load([basename '.cell_metrics.cellinfo.mat']);
       
        ntemp=0;
        for i = 1:numel(cell_metrics.putativeCellType) 
            if strcmp(cell_metrics.putativeCellType{i},'Pyramidal Cell') && ...
               strcmp(cell_metrics.brainRegion{i},'CA1')
               ntemp= ntemp+1;
            end
        end
        nCA1pyrSes(ses,1) = ntemp;
        nCA1pyr = nCA1pyr + ntemp;
        clear ntemp cell_metrics;
    end
end
