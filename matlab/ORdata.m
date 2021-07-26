%% Process OR dataset for SWR analysis 
clearvars;
dirData = 'A:\ORproject\';
animals = {'Wmaze2\OR15\','Wmaze2\OR18\','Wmaze3\OR22\','Wmaze3\OR21\','Wmaze3\OR23\'};% 'OR15\',  
sessions = {{'day1','day2','day3','day4','day10'},{'day1','day2','day3'},{'day1','day4','day3','day5'},...
            {'day2','day4'},{'day1','day5'}};

%% DONE
  % behavEpoch
  % cell_metrics
  % states
  % swr
  % tag cells
  % SWRpipeline     
        
% missing: brainRegions, pyr/int        
        
%% classify units
CXch = [29 35 25];
DGch = [47 16 46 19 45 18 44 20 42];
PFCch = 64:127;
for rat = 1:length(animals)
    for ses = 1:length(sessions{rat})
        cd([dirData animals{rat} sessions{rat}{ses}]);
        basename = bz_BasenameFromBasepath(pwd);   
        load([basename '.cell_metrics.cellinfo.mat']);
        
        for i = 1:numel(cell_metrics.UID)
            if ~strcmp('Pyramidal Cell',cell_metrics.putativeCellType{i}) && cell_metrics.firingRate(i)<3
               cell_metrics.putativeCellType{i}='Pyramidal Cell';
            end
        end
        
        for i = 1:numel(cell_metrics.UID)
            if find(CXch==cell_metrics.maxWaveformCh(i))
               cell_metrics.brainRegion{i} = 'NCX';
            elseif find(DGch==cell_metrics.maxWaveformCh(i))
               cell_metrics.brainRegion{i} = 'DGCA3'; 
            elseif find(PFCch==cell_metrics.maxWaveformCh(i))
               cell_metrics.brainRegion{i} = 'PFC';
            else
               cell_metrics.brainRegion{i} = 'CA1';
            end
        end
        
        save([basename '.cell_metrics.cellinfo.mat'],'cell_metrics');
        clear cell_metrics basename 
    end
end

%%
basename = bz_BasenameFromBasepath(pwd); 
load([basename '.session.mat']);
load([basename '.spikes.cellinfo.mat']);
spikes = loadSpikes('session',session,'forceReload',true);
cell_metrics = ProcessCellMetrics('session',session,'spikes',spikes,'excludeMetrics',{'monoSynaptic_connections','deepSuperficial'});
close all; clearvars; 

%%
for rat = 1:length(animals)
    for ses = 1:length(sessions{rat})
        cd([dirData animals{rat} sessions{rat}{ses}]);
        basename = bz_BasenameFromBasepath(pwd);   
        
        load([basename '.session.mat']);
        session.spikeSorting{1}.format = 'neurosuite';
        save([basename '.session.mat'],'session');
        
        load([basename '.spikes.cellinfo.mat']);
        spikes = loadSpikes('spikes',spikes,'session',session,'clusteringformat','neurosuite');
        cell_metrics = ProcessCellMetrics('session',session,'excludeMetrics',{'monoSynaptic_connections','deepSuperficial'});
        
        load([basename '.session.mat']);
        session.spikeSorting{1}.format = 'Phy';
        save([basename '.session.mat'],'session');   
        
        clear spikes basename cell_metrics session
        disp([animals{rat} sessions{rat}{ses} ' %%%%%%%%   DONE']);
    end
end
    
%%
for rat = 2:length(animals)
    for ses = 1:length(sessions{rat})
        try
        cd([dirData animals{rat} sessions{rat}{ses}]);
        basename = bz_BasenameFromBasepath(pwd); 
            
        session = sessionTemplate(pwd,'showGUI',false); %
        save([basename '.session.mat'],'session');
    
        %spikes = loadSpikes('session',session);
        cell_metrics = ProcessCellMetrics('session',session,'excludeMetrics',{'monoSynaptic_connections'});

        SleepScoreMaster(pwd,'noPrompts',true); % try to sleep score
        bz_thetaEpochs(pwd);           
        
        load('swrCh');
        ripples = bz_DetectSWR(swrCh+1,'saveMat',true);
        mkdir('Ripple_Profile');
        lfpRip = bz_GetLFP(swrCh(2));
        [wavAvg,lfpAvg]=bz_eventWavelet(lfpRip,ripples.peaks(1:500),'twin',[0.1 0.1]);
        saveas(gcf,['Ripple_Profile\ripWaveletSample.png']);
       
        clear cell_metrics session basename swrCh ripples lfpRip
        close all;
        disp([animals{rat} sessions{rat}{ses} ' %%%%%%%%   DONE']);
        catch
        warning([animals{rat} sessions{rat}{ses} ' %%%%%%%%   PROBLEM']);   
        end
   end
end
        