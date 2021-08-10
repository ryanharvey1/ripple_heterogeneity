% custom_for_girardeau

% has cell metrics but no .spk.
files = {'A:\Data\GirardeauG\Rat08\Rat08-20130708',...
    'A:\Data\GirardeauG\Rat08\Rat08-20130709'}

data_path = 'A:\Data\GirardeauG';
tic
files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat*']);
disp(toc)

data_path = 'A:\Data\GirardeauG';
tic
files = dir([data_path,'\**\*.spk.*']);
disp(toc)

sessions = unique({files.folder});

for i = 1:length(sessions)
    basepath = sessions{i};
    basename = bz_BasenameFromBasepath(basepath);
    disp(basepath)
    
    % check for needed files
    if exist(fullfile(basepath,[basename,'.spikes.cellinfo.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file') &&...
            exist(fullfile(basepath,[basename,'.SWRunitMetrics.mat']),'file')
        continue
    end
    
    % make sure there are the same n of res and spk
    if length(dir([basepath,filesep,'*.res.*'])) ~=...
            length(dir([basepath,filesep,'*.spk.*']))
        continue 
    end
    
    % check and make session.mat
    if ~exist([basename '.session.mat'],'file')
        session = sessionTemplate(basepath,'showGUI',false);
        session.epochs = get_kenji_epochs('basepath',basepath,'basename',basename);
        save(fullfile(basepath,[basename '.session.mat']),'session')
    else
        load(fullfile(basepath,[basename '.session.mat']))
    end

end