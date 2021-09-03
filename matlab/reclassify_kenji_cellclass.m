% reclassify_kenji_cellclass
data_path = 'A:\Data\Kenji';
tic
files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);
disp(toc)

for i = 1:length(files)
   basepath{i} = files(i).folder;
   basename{i} = basenameFromBasepath(files(i).folder);
end


load([data_path '\KenjiData2.mat']);
load([data_path '\ElePosition.mat']);

for s = 1:length(basepath)
    disp(basepath{s})
    load(fullfile(basepath{s},[basename{s},'.cell_metrics.cellinfo.mat']))
    
    for i = 1:size(ElePosition,1)
        if strcmp(ElePosition{i,2},basename{s})
            ses_code = ElePosition{i,5};
        end
    end

    % Get cell types
    cellSesMap=[];cleanSes=[];regionSes=[];iCellSes=cell(1,2); count=0;
    for i = 1:size(PyrIntMap.Map,1)
        if PyrIntMap.Map(i,18) == ses_code  
           count = count+1;
           cellSesMap(count,:) = PyrIntMap.Map(i,:);
           cleanSes(count,:) = Clean(i);
           iCellSes{1}(count,:) = iCell{1}(i); % 1=pyr
           iCellSes{2}(count,:) = iCell{2}(i);% 1= int
           regionSes(count,:) = Region(i);
        end    
    end
    putativeCellType = repmat({'unknown'},1,size(cell_metrics.UID,2));
    putativeCellType(iCellSes{1}==1) = {'Pyramidal Cell'};
    putativeCellType(iCellSes{2}==1) = {'Narrow Interneuron'};

    cell_metrics.putativeCellType = putativeCellType;
    
    save(fullfile(basepath{s},[basename{s},'.cell_metrics.cellinfo.mat']),'cell_metrics')
end



%% locate sessions
clear
data_path = 'A:\Data\Kenji';
tic
files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);
disp(toc)

for i = 1:length(files)
    basepath{i} = files(i).folder;
    basename{i} = basenameFromBasepath(files(i).folder);
%     load(fullfile(files(i).folder,files(i).name))
%     if size(cell_metrics.brainRegion,2) == 0
%         disp(files(i).folder)
%         cell_metrics.brainRegion = repmat({'unknown'},1,size(cell_metrics.UID,2));
%         save(fullfile(files(i).folder,files(i).name),'cell_metrics')
%     end
end


%% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepath,'basenames',basename);

%% classify shank region

df = readtable('A:\Data\Kenji\ElePosition.Mizuseki_v2.csv');
sessions = unique(cell_metrics.sessionName(strcmp(cell_metrics.brainRegion,'unknown')));
for i = 1:length(sessions)
    cell_metric_idx = strcmp(cell_metrics.sessionName,sessions{i});

    brain_region = table2cell(df(strcmp(df.session,sessions{i}),...
        cell_metrics.shankID(cell_metric_idx)+4));
    
    cell_metrics.brainRegion(cell_metric_idx) = brain_region;
end
cell_metrics.brainRegion = upper(cell_metrics.brainRegion);

%% make custom naming for kenji dataset
ids = {'ec013','ec014','ec016','f01_m',...
      'g01_m','gor01','i01_m','j01_m','km01',...
      'nlx'};

% for id = ids
%     idx = contains(cell_metrics.sessionName,id);
%     cell_metrics.animal(idx) = id;
% end
% need to save manually... :( cell explorer is not saving the new names
for i = 1:length(basepath)
    disp(basepath{i})
    load(fullfile(basepath{i},[basename{i},'.cell_metrics.cellinfo.mat']))
    for id = ids
        idx = contains(cell_metrics.sessionName,id);
        cell_metrics.animal(idx) = id;
    end
    save(fullfile(basepath{i},[basename{i},'.cell_metrics.cellinfo.mat']),'cell_metrics')
end

%%
cell_metrics = CellExplorer('metrics',cell_metrics);

% basename = basenameFromBasepath(basepath);

%%
X = cell_metrics.waveforms.filt_zscored';
y = cell_metrics.putativeCellType';
X = X(~strcmp(y,'unknown'),:);
y = y(~strcmp(y,'unknown'));

Mdl = fitclinear(X,y);

X = cell_metrics.waveforms.filt_zscored';
y = cell_metrics.putativeCellType';
X = X(strcmp(y,'unknown'),:);
labels = predict(Mdl,X);

waveforms = cell_metrics.waveforms.filt_zscored(:,strcmp(y,'unknown'));
figure
plot(waveforms(:,contains(labels,'Pyramidal Cell')))
figure
plot(waveforms(:,contains(labels,'Narrow Interneuron')))


cell_metrics.putativeCellType(strcmp(y,'unknown')) = labels;

%% make deep sup classes
clear
data_path = 'A:\Data\Kenji';
tic
files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);
disp(toc)

for i = 1:length(files)
    basepath{i} = files(i).folder;
    basename{i} = basenameFromBasepath(files(i).folder);
end
for i = 1:length(basepath)
    disp(basepath{i})

    if ~exist(fullfile(basepath{i},[basename{i},'.deepSuperficialfromRipple.channelinfo.mat']))
%             load(fullfile(basepath{i},[basename{i},'.spikes.cellinfo.mat']))
    load(fullfile(basepath{i},[basename{i},'.session.mat']))
%     load(fullfile(basepath{i},[basename{i},'.cell_metrics.cellinfo.mat']))
        disp('running classification_DeepSuperficial...')
        deepSuperficialfromRipple = classification_DeepSuperficial(session);
    end
%     cell_metrics = ProcessCellMetrics('basepath',basepath{i},...
%         'showGUI',false,...
%         'spikes',spikes,...
%         'getWaveformsFromDat',false,...
%         'manualAdjustMonoSyn',false,...
%         'session',session);
end