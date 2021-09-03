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

%%


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
