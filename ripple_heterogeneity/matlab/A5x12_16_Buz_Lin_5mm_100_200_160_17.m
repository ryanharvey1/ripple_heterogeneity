
%% % % % % % % % % % % % % % % % % % % % % % % 
% A5x12-16-Buz-Lin-5mm-100-200-160-177 (64 ch, 5 shanks, poly 2  Custom)
% % % % % % % % % % % % % % % % % % % % % % 
clear chanCoords
chanCoords.probe = 'A5x12-16-Buz-Lin-5mm-100-200-160-177';
chanCoords.source = 'Manually entered';
chanCoords.layout = 'poly 2';
chanCoords.nShanks = 5;
chanCoords.nChannels = 64;
chanCoords.shankSpacing = 160; % in um
chanCoords.verticalSpacing = 20; % in um
chanCoords.horizontalSpacing = 17.32; % in um

coords = readtable("D:\github\ripple_heterogeneity\data\electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177 (1).csv");

% % % % % % % % % % % % % % % % % % % % % % 
% Generating map
chanCoords.x = coords.x_um_;
chanCoords.y = coords.y_um_;

shankSpacing = chanCoords.shankSpacing;
verticalSpacing = chanCoords.verticalSpacing;
horizontalSpacing = chanCoords.horizontalSpacing;

% % % % % % % % % % % % % % % % % % % % % % 
% Adding processing info
chanCoords.processinginfo.function = 'ProcessCellMetrics';
chanCoords.processinginfo.date = now;
try
    chanCoords.processinginfo.username = char(java.lang.System.getProperty('user.name'));
    chanCoords.processinginfo.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
catch
    disp('Failed to retrieve system info.')
end

%% % % % % % % % % % % % % % % % % % % % % % 
% Plotting and saving channel coordinates
% % % % % % % % % % % % % % % % % % % % % % 
x_range = range(chanCoords.x);
y_range = range(chanCoords.y);
if x_range > y_range
    fig_width = 1600;
    fig_height = ceil(fig_width*y_range/x_range)+200;
else
    fig_height = 1000;
    fig_width = ceil(fig_height*x_range/y_range)+200;
end
fig1 = figure('Name',['Channel map: ' chanCoords.probe ],'position',[5,5,fig_width,fig_height]); movegui(fig1,'center')
plot(chanCoords.x,chanCoords.y,'.k'), hold on
text(chanCoords.x,chanCoords.y,num2str([1:numel(chanCoords.y)]'),'VerticalAlignment', 'bottom','HorizontalAlignment','center');
title({' ','Channel map',' '}), xlabel('X (um)'), ylabel('Y (um)'),
