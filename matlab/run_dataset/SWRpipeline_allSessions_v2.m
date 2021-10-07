% SWRpipeline_allSessions_v2

animal = {'AB1','AB3','AB4','AYA4','AYA6','AYA7','AYA9','AYA10',...
    'OML5','OML3','OML7','OML8','OML10','OML18','OML19',...
    'Wmaze2\OR15','Wmaze2\OR18','Wmaze3\OR22','Wmaze3\OR21','Wmaze3\OR23',...
    'GrosmarkAD\Cicero','GrosmarkAD\Buddy','GrosmarkAD\Achilles','GrosmarkAD\Gatsby',...
    'Kenji'};

dataDir1 = 'A:\Data\';
dataDir2 = 'A:\OptoMECLEC\';
dataDir3 = 'A:\ORproject\';

for a = length(animal)
    if strncmp('OML',animal{a},3)
        base_path = dataDir2;
    elseif strncmp('Wmaze',animal{a},5)
        base_path = dataDir3;
    else
        base_path = dataDir1;
    end
    files = dir([base_path,...
        animal{a},...
        filesep,'**',filesep,...
        filesep,'**',filesep,...
        '*.cell_metrics.cellinfo.mat']);
    for f = 1:length(files)
        basepath = files(f).folder;
        basename = bz_BasenameFromBasepath(basepath);
        if exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file') &&...
                exist(fullfile(basepath,[basename,'.SWRunitMetrics.mat']),'file')
            load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
            load(fullfile(basepath,[basename,'.SWRunitMetrics.mat']))

        end
                
    end
end

