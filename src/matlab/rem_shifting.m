% rem_shifting

% These are the criteria for defining REM shifting cell from Mizuseki, 2011. 
% Neurons with <120° or >300° preferred theta phases during REM (magenta)
% were designated as REM-shifting cells, whereas those between 120° to 300°
% (blue) were designated as nonshifting cells.

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

date_thres = '8-Dec-2021';
basepaths = unique(df.basepath);
for i = 1:length(basepaths)
    % make input
%     get_rem_shift_input(basepaths{i})
    disp(basepaths{i})
    basename = basenameFromBasepath(basepaths{i});
    
    if exist(fullfile(basepaths{i},[basename,'.theta_rem_shift.mat']),'file')
        % check date
        file = dir(fullfile(basepaths{i},[basename,'.theta_rem_shift.mat']));
        if file.datenum >= datenum(date_thres)
            continue
        end
        % remove file first
        delete(fullfile(basepaths{i},[basename,'.theta_rem_shift.mat']))
    end
    get_rem_shift('basepath',basepaths{i});
end





% 
% basepath = unique(df.basepath(contains(df.basepath,'Kenji')));
% basepath = basepath{9}
% % basepath = df.basepath{1};
% basename = basenameFromBasepath(basepath);
% 
% load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))

% % set up for finding theta channel for kenji
% ca1_shanks = unique(cell_metrics.shankID(contains(lower(cell_metrics.brainRegion),'ca1')));
% ca1_channels = [cell_metrics.general.electrodeGroups{ca1_shanks}];
% 
% theta_channel = load(fullfile(basepath,[basename,'.Theta_eegCh']));
% 
% theta_channel = intersect(ca1_channels,theta_channel(:,2));

%% set up for finding theta channel for Girardeau
% files = dir(fullfile(basepath,'*PhaseModulation*.mat'));
% files = {files.name};
% files = files(contains(files,'Hpc') & contains(files,'Theta'));
% data = load(fullfile(basepath,files{1}));
% theta_channel = data.info.channel;
%%
% spikes = loadSpikes('basepath',basepath,'basename',basename);
% [lfp] = getLFP(theta_channel,'basepath',basepath);
% 
% load(fullfile(basepath,[basename,'.SleepState.states.mat']));
%%
% figure;
% 
% states = fields(SleepState.ints);
% viewwin = [SleepState.detectorinfo.StatePlotMaterials.t_clus(1),...
%     SleepState.detectorinfo.StatePlotMaterials.t_clus(end)];
% colors = lines(length(states));
% for i = 1:length(states)
%     disp(states{i})
%     plot(SleepState.ints.(states{i})',(-i)*ones(size(SleepState.ints.(states{i})))',...
%         'color',colors(i,:),'LineWidth',8)
%     hold on;
% end
% xlim([viewwin]) 
% ylim([-i-1 0])
% set(gca,'YTick',[-i:-1])
% set(gca,'YTickLabel',flipud(states))

%%
% wake_n = find(strcmp(SleepState.idx.statenames,'WAKE'));
% idx_wake = SleepState.idx.states == wake_n;
% 
% rem_n = find(strcmp(SleepState.idx.statenames,'REM'));
% idx_rem = SleepState.idx.states == rem_n;
% 
% theta_n = find(strcmp(SleepState.idx.theta_epochs.statenames,'THETA'));
% idx_theta = SleepState.idx.theta_epochs.states == theta_n;
%   
% idx_wake_theta = idx_wake & idx_theta;
% idx_rem_theta = idx_theta & idx_rem;

    
% [PhaseLockingData_wake] = phaseModulation(spikes,...
%                                     lfp,...
%                                     [6 12],...
%                                     'saveMat', false,...
%                                     'plotting',false,...
%                                     'numBins',60,...
%                                     'basepath',basepath,...
%                                     'intervals',ToIntervals(idx_wake),...
%                                     'powerThresh',0);
%                                 
% [PhaseLockingData_rem] = phaseModulation(spikes,...
%                                     lfp,...
%                                     [6 12],...
%                                     'saveMat', false,...
%                                     'plotting',false,...
%                                     'numBins',60,...
%                                     'basepath',basepath,...
%                                     'intervals',ToIntervals(idx_rem),...
%                                     'powerThresh',0);
                                
% Neurons with <120° or >300° preferred theta phases during REM (magenta)
% were designated as REM-shifting cells, whereas those between 120° to 300°

% angles = rad2deg(PhaseLockingData_rem.phasestats.m);
% rem_shift = [angles<120 | angles>300] & PhaseLockingData_rem.phasestats.p<0.05;
% non_rem_shift = [angles>120 & angles<300] & PhaseLockingData_rem.phasestats.p<0.05;
% 
% figure;
% plot(PhaseLockingData_wake.phasebins,...
%     PhaseLockingData_wake.phasedistros(:,rem_shift))
% 
% figure;
% plot(PhaseLockingData_wake.phasebins,...
%     PhaseLockingData_wake.phasedistros(:,non_rem_shift))

%%
% figure;                            
% plot(PhaseLockingData_wake.phasebins,...
%     PhaseLockingData_wake.phasedistros(:,PhaseLockingData_wake.phasestats.p<0.05))
% 
% 
% figure;                            
% plot(PhaseLockingData_rem.phasebins,...
%     PhaseLockingData_rem.phasedistros(:,PhaseLockingData_rem.phasestats.p<0.05))

