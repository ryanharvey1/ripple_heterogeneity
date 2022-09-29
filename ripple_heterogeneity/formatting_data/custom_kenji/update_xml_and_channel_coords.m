
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

df = df(contains(df.basepath,'ec016'),:)
xml_file = 'Z:\Data\Kenji\ec016.100_121\ec016.100_121.xml'

for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    if exist(fullfile(basepath,[basename,'_old_old.xml']),'file')
        continue
    end
    disp(basepath)
    movefile(fullfile(basepath,[basename,'.xml']), fullfile(basepath,[basename,'_old_old.xml']))
    
    copyfile(xml_file, fullfile(basepath,[basename,'.xml']));
    
    load(fullfile(basepath,[basename,'.session.mat']),'session')
    
    session = import_xml2session(xml_file,session);
    
    session.extracellular.chanCoords.layout = 'staggered';
    session.extracellular.chanCoords.shankSpacing = 200;
    session.extracellular.chanCoords.verticalSpacing = 20;
    
    if isfield(session.animal,'probeImplants')
        session.animal.probeImplants{1}.verticalSpacing = 20;
    end
    chanCoords = generateChanCoords(session);
    session.extracellular.chanCoords = chanCoords;
    
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end

%%
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

df = df(contains(df.basepath,'ec016'),:);

for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    disp(basepath)
    
    load(fullfile(basepath,[basename,'.session.mat']),'session')
    
    if isfield(session.extracellular,'chanCoords')
        if ~isempty(session.extracellular.chanCoords.x)
            continue
        end
    end
    
    session.extracellular.chanCoords.layout = 'staggered';
    session.extracellular.chanCoords.shankSpacing = 200;
    session.extracellular.chanCoords.verticalSpacing = 20;
    
    chanCoords = generateChanCoords(session);
    session.extracellular.chanCoords = chanCoords;
    
    save(fullfile(basepath,[basename,'.session.mat']),'session')
    
end

%%
%     function plotChannelMap1(~,~)
%         readBackChanCoords
%         if isfield(session,'extracellular') && isfield(session.extracellular,'chanCoords')
%             chanCoords = session.extracellular.chanCoords;
%             x_range = range(chanCoords.x);
%             y_range = range(chanCoords.y);
%             if x_range > y_range
%                 fig_width = 1600;
%                 fig_height = ceil(fig_width*y_range/x_range)+200;
%             else
%                 fig_height = 1000;
%                 fig_width = ceil(fig_height*x_range/y_range)+200;
%             end
%             fig1 = figure('Name','Channel coordinates','position',[5,5,fig_width,fig_height],'visible','off'); movegui(fig1,'center')
%             ax1 = axes(fig1);
%             plot(ax1,chanCoords.x,chanCoords.y,'.k'), hold on
%             text(ax1,chanCoords.x,chanCoords.y,num2str([1:numel(chanCoords.x)]'),'VerticalAlignment', 'bottom','HorizontalAlignment','center');
%             title(ax1,{' ','Channel coordinates',' '}), xlabel(ax1,'X (um)'), ylabel(ax1,'Y (um)')
%             set(fig1,'visible','on')
%         else
%             MsgLog('No channel coords data available',4)
%         end
%     end