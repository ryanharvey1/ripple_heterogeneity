% fix_cell_explorer_var_dims

animal = {'AB1','AB3','AB4','AYA4','AYA6','AYA7','AYA9','AYA10',...
    'OML5','OML3','OML7','OML8','OML10','OML18','OML19',...
    'Wmaze2\OR15','Wmaze2\OR18','Wmaze3\OR22','Wmaze3\OR21','Wmaze3\OR23',...
    'GrosmarkAD\Cicero','GrosmarkAD\Buddy','GrosmarkAD\Achilles','GrosmarkAD\Gatsby',...
    'Kenji'};

dataDir1 = 'A:\Data\';
dataDir2 = 'A:\OptoMECLEC\';
dataDir3 = 'A:\ORproject\';

for a = 1:length(animal)
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
        if exist(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'file')
            load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
            
            flag = [];
            for field = fields(cell_metrics)'
                [r,c] = size(cell_metrics.(field{1}));
                if r>c
                    flag = [flag;field];
                end
            end
            for field = flag
                cell_metrics.(field{1}) = cell_metrics.(field{1})';
            end
            if ~isempty(flag)
                save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
            end
        end
    end
end