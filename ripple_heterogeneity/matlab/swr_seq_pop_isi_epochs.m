
% basepath = 'A:\Data\GrosmarkAD\Cicero\Cicero_09172014';

savepath = 'D:\projects\ripple_heterogeneity\swr_seq_isi_low_high_partic';

df = readtable('D:\projects\ripple_heterogeneity\sessions.csv');
df_meta = readtable('D:\projects\ripple_heterogeneity\session_metadata.csv');

idx = strcmp(df_meta.opto,'none') &...
    [strcmp(df_meta.task,'t_maze') |...
    strcmp(df_meta.task,'cheesboard') |...
    strcmp(df_meta.novel,'novel')];
df_meta = df_meta(idx,:);


for i = 1:length(df_meta.animal)
    idx_row(i) = find(contains(df.basepath,df_meta.animal{i}) & contains(df.basename,df_meta.day{i}));
end

basepaths = unique(df.basepath(idx_row,:));

% basepaths = basepaths(contains(basepaths,'GrosmarkAD'));

particip_thres = [0.096, 0.25]; % based on pre task grosmark pyr ca1

epochs = {'pre','task','post'};
for ep = 1:length(epochs)
    savepath_ = fullfile(savepath,epochs{ep});
    if ~exist(savepath_,'dir')
        mkdir(savepath_)
    end
    for basepath = basepaths'
        basename = basenameFromBasepath(basepath{1});
        load([basepath{1},filesep,basename,'.session.mat'])
        load(fullfile(basepath{1},[basename '.SWRunitMetrics.mat']));
        
        if length(session.epochs) ~= 3
            continue
        end
        
        custom_grouping = repmat({'unknown'},1,length(SWRunitMetrics.(epochs{ep}).particip));
        custom_grouping(SWRunitMetrics.(epochs{ep}).particip < particip_thres(1)) = {'low_particip'};
        custom_grouping(SWRunitMetrics.(epochs{ep}).particip > particip_thres(2)) = {'high_particip'};

        restrict_ts = [session.epochs{ep}.startTime,session.epochs{ep}.stopTime];

        swr_seq_pop_isi('basepath',basepath{1},'savepath',savepath_,...
            'parallel',false,'restrict_ts',restrict_ts,...
            'custom_grouping',custom_grouping,...
            'binary_class_variable','particip',...
            'grouping_names',{'low_particip','high_particip'})
    end
end

% restrict_ts = [session.epochs{1, 1}.startTime,session.epochs{1, 1}.stopTime];

% 'AYA9'	'day20'	'none'
% 'AYA10'	'day25'	'none'
% 'AYA10'	'day27'	'none'
% 'AYA10'	'day31'	'none'
% 'AYA10'	'day32'	'none'
% 'AYA10'	'day34'	'none'

% epochs = {'PREEpoch','MazeEpoch','POSTEpoch'};
% for basepath = basepaths'
%     load([basepath{1},filesep,basenameFromBasepath(basepath{1}),'.session.mat'])
%     load([basepath{1},filesep,basenameFromBasepath(basepath{1}),'_sessInfo.mat'])
%     for e = 1:length(session.epochs)
%         if ~isfield(session.epochs{e},'startTime')
%             session.epochs{e}.startTime = sessInfo.Epochs.(epochs{e})(1);
%             session.epochs{e}.stopTime = sessInfo.Epochs.(epochs{e})(2);   
%         end
%     end
%     save([basepath{1},filesep,basenameFromBasepath(basepath{1}),'.session.mat'],'session')
% end
