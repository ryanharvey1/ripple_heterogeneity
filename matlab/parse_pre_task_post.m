function parse_pre_task_post(session,basepath,basename,ripples,spikes)

for i = 1:length(session.epochs)
    epoch{i} = session.epochs{i}.name;
end
disp(epoch')
% idx = find(diff(find(contains(epoch,'sleep')))>1)
% pre = session.epochs{i}

epochs = {'pre','task','post'};

% x = input('pre/task/post ? (y/n)','s');
% if strcmp(x,'n')
%     x = input('list epoch','s');
% end
% if exist(fullfile(basepath,[basename '.SWRunitMetrics.mat']),'file')
%     
% end

epoch_struct = [];
[pre_idx,task_idx,post_idx ] = locate_idx_pre_task_post(session);
if ~isempty(pre_idx)
    epoch_struct.pre = [session.epochs{pre_idx(1)}.startTime,...
        session.epochs{pre_idx(2)}.stopTime];
end
if ~isempty(task_idx)
    epoch_struct.task = [session.epochs{task_idx(1)}.startTime,...
        session.epochs{task_idx(2)}.stopTime];
end
if ~isempty(post_idx)
    epoch_struct.post = [session.epochs{post_idx(1)}.startTime,...
        session.epochs{post_idx(2)}.stopTime];
end
% disp('manual for now...sorry')
% disp('manipulate & run epoch_struct and type dbcont')
% keyboard
% if isempty(epoch_struct)
%     epoch_struct.pre = [session.epochs{1}.startTime,session.epochs{2}.stopTime];
%     epoch_struct.task = [session.epochs{3}.startTime,session.epochs{5}.stopTime];
%     epoch_struct.post = [session.epochs{6}.startTime,session.epochs{7}.stopTime];
%     
%     epoch_struct.task_2 = [session.epochs{8}.startTime,session.epochs{14}.stopTime];
%     epoch_struct.post_2 = [session.epochs{15}.startTime,session.epochs{17}.stopTime];
%     
%     epoch_struct.task_3 = [session.epochs{17}.startTime,session.epochs{19}.stopTime];
%     epoch_struct.post_3 = [session.epochs{20}.startTime,session.epochs{21}.stopTime];
% end
SWRunitMetrics = [];
SWRunitMetrics.epoch_struct = epoch_struct;

if isfield(epoch_struct,'pre')
    SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
        basename,...
        ripples,...
        spikes,...
        SWRunitMetrics,...
        epochs,...
        epoch_struct,...
        'pre',...
        'pre');
end

if isfield(epoch_struct,'task')
    SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
        basename,...
        ripples,...
        spikes,...
        SWRunitMetrics,...
        epochs,...
        epoch_struct,...
        'task',...
        'task');
end
if isfield(epoch_struct,'post')
    SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
        basename,...
        ripples,...
        spikes,...
        SWRunitMetrics,...
        epochs,...
        epoch_struct,...
        'post',...
        'post');
end
save(fullfile(basepath,[basename '.SWRunitMetrics.mat']),'SWRunitMetrics')
end

function SWRunitMetrics = calc_sliced_unitSWRmetrics(basepath,...
                                                    basename,...
                                                    ripples,...
                                                    spikes,...
                                                    SWRunitMetrics,...
                                                    epochs,...
                                                    epoch_struct,...
                                                    epoch_search,...
                                                    epoch_label)
                                                
start_end = epoch_struct.(epochs{contains(epochs,epoch_search)});
if isfield(ripples,'times')
    idx = ripples.times(:,1) >= start_end(1) &...
        ripples.times(:,2) <= start_end(2);
    
    ripSpk = bz_getRipSpikes('basepath',basepath,...
                            'basename',basename,...
                            'spikes',spikes,...
                            'events',ripples.times(idx,:),...
                            'saveMat',false); 
                        
elseif isfield(ripples,'timestamps')
    idx = ripples.timestamps(:,1) >= start_end(1) &...
        ripples.timestamps(:,2) <= start_end(2);
    
    ripSpk = bz_getRipSpikes('basepath',basepath,...
                            'basename',basename,...
                            'spikes',spikes,...
                            'events',ripples.timestamps(idx,:),...
                            'saveMat',false); 
end

SWRunitMetrics.(epoch_label) = bz_unitSWRmetrics(ripSpk);

end

function [pre_idx,task_idx,post_idx] = locate_idx_pre_task_post(session)
% locate the first pre/task/post epochs. 
% They can be contiguous as in (sleep,sleep,task,sleep,sleep), and 
% the below will figure out where the epochs should go

for i = 1:length(session.epochs)
    epoch{i,1} = lower(session.epochs{i}.name);
end
for i = 1:length(epoch)
   epoch_temp = strsplit(epoch{i},'_');
   epoch{i,1} = epoch_temp{2};
end

idx = contains(epoch,'sleep');
clear epoch_label
for i = 1:length(epoch)
    if i==1
        if contains(epoch{i},'sleep') 
            epoch_label{i,1} = 'pre';
        else
            epoch_label{i,1} = 'task';
        end
    else
        if contains(epoch{i},'sleep') &&...
                ~contains(epoch_label{i-1},'task') &&...
                ~contains(epoch_label{i-1},'post')
            epoch_label{i,1} = 'pre';
        elseif contains(epoch{i},'sleep') &&...
                [contains(epoch_label{i-1},'task') ||...
                contains(epoch_label{i-1},'post')]
            epoch_label{i,1} = 'post';
        else
            epoch_label{i,1} = 'task';
        end
    end
end

pre_idx = [];
task_idx = [];
post_idx = [];

[start,ends,ngroups]=find_groups(contains(epoch_label,'pre'));
if ngroups>0
    pre_idx = [start(1),ends(1)];
end
[start,ends,ngroups]=find_groups(contains(epoch_label,'task'));
if ngroups>0
    task_idx = [start(1),ends(1)];
end
[start,ends,ngroups]=find_groups(contains(epoch_label,'post'));
if ngroups>0
    post_idx = [start(1),ends(1)];
end
end