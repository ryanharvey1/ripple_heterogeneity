

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
% df.animal = repmat({'unknown'},length(df.basepath),1);
% for i = 1:length(df.basepath)
%     basepath = df.basepath{i};
%     basename = basenameFromBasepath(basepath);
%     load(fullfile(basepath,[basename,'.session.mat']))
%     df.animal(i) = {session.animal.name};
% end
% unique(df.animal)
%%

basepaths = df.basepath(contains(df.basepath,'Can'));
% load manually fixed session
% this_session = load('day03.session.mat')
basepath = basepaths{1};
basename = basenameFromBasepath(basepath);
this_session = load(fullfile(basepath,[basename,'.session.mat']))

% iter over other sessions for this animal
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.session.mat']))
    
    session.animal.probeImplants = this_session.session.animal.probeImplants;
    session.extracellular.chanCoords = this_session.session.extracellular.chanCoords;
    session.analysisTags = this_session.session.analysisTags;
    
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end


%%

session.extracellular.chanCoords.y .* deepSuperficial_ChDistance3'


for jj = 1:session.extracellular.nElectrodeGroups
    % Analysing the electrode groups separately
%     fprintf(['Analysing electrode group ', num2str(jj),', ']);
    
    % Get list of channels belong to electrode group (1-indexed)
    ripple_channels{jj} = session.extracellular.electrodeGroups.channels{jj};
    
    % remove ripple channels that are labeled 'Bad' 
    ripple_channels{jj}(ismember(ripple_channels{jj},channels_to_exclude)) = [];
%             disp(deepSuperficial_ChDistance(ripple_channels{jj}))

    if any(deepSuperficial_ChDistance(ripple_channels{jj}) > 0) & any(deepSuperficial_ChDistance(ripple_channels{jj}) < 0)
        disp(deepSuperficial_ChDistance(ripple_channels{jj}))
    end
end

x = session.extracellular.chanCoords.x(ripple_channels{jj})
y = session.extracellular.chanCoords.y(ripple_channels{jj})

scatter(x,y)


ia2_test = abs(session.extracellular.chanCoords.y(ripple_channels{jj}))

interp1(SWR_diff2([indx,indx+1]),[ia2_test(indx),ia2_test(indx+1)],0)

unique(cell_metrics.deepSuperficial(isnan(cell_metrics.deepSuperficialDistance)))