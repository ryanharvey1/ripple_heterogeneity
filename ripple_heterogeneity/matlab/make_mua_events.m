
df = readtable("Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv");

parfor i = 1:length(df.basepath)
    disp(df.basepath{i})
    
    run(df.basepath{i})
end

function run(basepath)
    basename = basenameFromBasepath(basepath);
    save_file = fullfile(basepath,[basename,'.mua_ca1_pyr.events','.mat']);
    if exist(save_file,'file')
       return 
    end
    spikes = importSpikes('basepath',basepath,...
        'brainRegion',"CA1",...
        'cellType',"Pyr");
    if isempty(spikes.UID)
        return
    end
    find_HSE(spikes,'basepath',basepath,...
        'name','mua_ca1_pyr',...
        'save_classic_event_log',false);
end