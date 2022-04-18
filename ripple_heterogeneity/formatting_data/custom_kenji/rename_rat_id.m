% rename_rat_id
% make custom naming for kenji dataset
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'Kenji'),:);

ids = {'2006-4-10','2006-4-18','2006-6-7','2006-6-12','2006-6-13',...
    'ec013','ec014','ec016','f01_m',...
    'g01_m','gor01','i01_m','j01_m','km01',...
    'nlx'};

for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath,[basename,'.session.mat']))
    
    for id = 1:length(ids)
        name = ids{id};
        if contains(cell_metrics.sessionName{1},name)
            cell_metrics.animal = repmat({name},1,length(cell_metrics.animal));
            
            session.animal.name = name;
            
            save(fullfile(basepath,[basename,'.session.mat']),'session')
            save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
            break
        end
    end 
end