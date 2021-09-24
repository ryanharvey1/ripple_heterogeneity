% add_participation_to_cell_metrics

df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepaths = unique(df.basepath);

% find sessions to run
% run_it_idx = zeros(length(basepaths),1);
% parfor i = 1:length(basepaths)
%     basepath = basepaths{i};
%     basename = basenameFromBasepath(basepath);
%     run_it_idx(i) = check_file(basepath,basename);
% end
% basepaths = basepaths(logical(run_it_idx));

basepaths = flipud(basepaths);

WaitMessage = parfor_wait(length(basepaths));
parfor i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    main(basepath,basename)
    WaitMessage.Send;
end
WaitMessage.Destroy;


function run_it = check_file(basepath,basename)
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
run_it = ~isfield(cell_metrics,'ripple_particip');
end

function main(basepath,basename)
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))

if isfield(cell_metrics,'ripple_particip')
    return
end
load(fullfile(basepath,[basename,'.session.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
parameters = cell_metrics.general.processinginfo.params;

cell_metrics = customCalculations.('swr_unit_metrics')(cell_metrics,session,{spikes},parameters);

save_new_metrics(basepath,basename,cell_metrics)

end
function save_new_metrics(basepath,basename,cell_metrics)
dirname = 'revisions_cell_metrics';
saveAs = 'cell_metrics';
if ~(exist(fullfile(basepath,dirname),'dir'))
    mkdir(fullfile(basepath,dirname));
end
if exist(fullfile(basepath,[basename,'.',saveAs,'.cellinfo.mat']),'file')
    copyfile(fullfile(basepath,[basename,'.',saveAs,'.cellinfo.mat']),...
        fullfile(basepath, dirname, [saveAs, '_',datestr(clock,'yyyy-mm-dd_HHMMSS'), '.mat']));
end
save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end