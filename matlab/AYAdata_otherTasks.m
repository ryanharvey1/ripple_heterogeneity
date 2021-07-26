%% 
% this sessions are with different task than linear track (cheesboard, t maze, etc.) 

dataDir = 'A:\Data\';
animal = {'AB3','AB4','AYA7','AYA9','AYA10'};
day = {{'AB3_47_49','AB3_50_51','AB3_55_57','AB3_60'},{'day03','day07','day08','day09','day11'},...
       {'day22','day24','day25','day30'},{'day15','day16','day20'},{'day25','day27','day31','day32','day34'}};

%%
for a = 2:numel(animal)
    for d = 5:numel(day{a})
        cd([dataDir animal{a} '\' day{a}{d}]);
        basename = bz_BasenameFromBasepath(pwd);
        
        SleepScoreMaster(pwd,'noPrompts',true); % try to sleep score
        bz_thetaEpochs(pwd);          
    end
end

% Problems with state score:
    % 'AB3_50_51' 'AB3_47_49'
    
% Repeat cell_metrics (waveCh1 es 0)            
    
%%





