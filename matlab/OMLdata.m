%% Process OML dataset for SWR analysis 

dirData = 'A:\OptoMECLEC\';

animals = {'OML5','OML3','OML7','OML8','OML10'};
sessions = {{'day4','day5','day6','day8','day9','day12','day13','day20','day21'},...
            {'day11','day12','day13','day14','day16','day17'},...
            {'day6','day7','day9'},...
            {'day4','day5','day6','day7','day8','day9','day17','day20'},...
            {'day5','day7','day9'}}; 
condition ={[1 0 1 1 0 8 3 0 1],...
               [1 0 1 1 3 3],...
               [1 0 1],...
               [0 2 0 2 0 2 0 0 2 0],...
               [2 2 2]}; 
% 0=sham, 1=MEC silencing; 2=LEC silencing; 3=prb stm; 8=mierda

%% DONE
  % behavEpoch
  % cell_metrics
  % states
  % swr
  % tag cells
  % SWRpipeline
  
for rat = 1:length(animals)
    for ses = 1:length(sessions{rat})
        try
            basepath = [dirData animals{rat} '\' sessions{rat}{ses}];
            cd(basepath);
            basename = bz_BasenameFromBasepath(pwd);
            load([basename '.spikes.cellinfo.mat']);
            load([basename '.session.mat']);
            
            if exist([basename '.cell_metrics.cellinfo.mat'],'file')
               load([basename '.cell_metrics.cellinfo.mat']);
            else
               cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
            end
            
            load([basename '.sessionInfo.mat']);
            cell_metrics.brainRegion = spikes.region;
            
            for i = 1:numel(cell_metrics.UID)
                if ~strcmp('Pyramidal Cell',cell_metrics.putativeCellType{i}) && cell_metrics.firingRate(i)<5
                   cell_metrics.putativeCellType{i}='Pyramidal Cell';
                end
            end
            
            save([basename '.cell_metrics.cellinfo.mat'],'cell_metrics');
            
            close all
            clear cell_metrics spikes session sessionInfo
            
            disp(['%%%%% ' animals{rat} sessions{rat}{ses} '  DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
        catch
            warning(['%%%%% ' animals{rat} sessions{rat}{ses} '  ERROR %%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
            
        end
    end
end


