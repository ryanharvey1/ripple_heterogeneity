
% add_familiarity_to_session
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'Kenji'),:);

load('Z:\Data\Kenji\KenjiData3.mat')

for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    
    disp(basepath)
    
    load(fullfile(basepath,[basename,'.session.mat']))
    
    temp_beh = Beh(contains(Beh(:,2),basename),:);
    
    for ep = 1:length(session.epochs)
       
        if ~any(strfind(session.epochs{ep}.name,temp_beh{ep,5}))
           disp('could not find session indicator') 
        end
        
        session.epochs{ep}.environment = temp_beh{ep,5};
        session.epochs{ep}.behavioralParadigm = str2double(temp_beh{ep,6});
    end
    
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end

