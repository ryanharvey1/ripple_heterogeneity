

figure;scatter(chanCoords.x,chanCoords.y,[],1:length(chanCoords.y))

% visualize channels with chan numbers
% figure;scatter(coords.x_um_,coords.y_um_,[],1:length(coords.y_um_))
hold on
text(chanCoords.x,chanCoords.y,string(1:length(chanCoords.x)),...
    'VerticalAlignment','bottom','HorizontalAlignment','center','fontsize',8);


five_by_12_y = [-100, -80, -60, -40, -20, 0, -120, -140, -150, -130, -10, -30,...
    -50, -70, -90, -110, -20, -40, -120, 0, -140, -60, -130, -150, -10, -80,...
    -50, -30, -70, -100, -110, -90, -100, -80, -60, -110, -40, -20, 0, -90,...
    -120, -140, -150, -70, -130, -10, -30, -50, -80, -100, -40, -60, 0, -20,...
    -140, -120, -130, -150, -30, -10, -70, -50, -110, -90, -80, -100, -110,...
    -10, -90, -10, -30, -30, -50, -70, -100, -140, -80, -60, -40, -120, -40,...
    -60, 0, -20, -100, -80, -90, -110, -30, -10, -70, -50, -40, -60, 0, -20,...
    -10, -30, -50, -70, -60, -40, -20, 0, -80, -100, -110, -90, -10, -30, -50,...
    -70, -50, -130, -90, -70, -110, -150, -40, -60, -20, -20, -80, 0, -100, 0,...
    -90, -110, -10, -30, -50, -70, -60, -40, -20, 0, -80, -100, -110, -90, -10,...
    -30, -50, -70, -80, -60, -120, -100, -140, -40, -40, -60, -20, -20, -80, 0,...
    -100, 0, -90, -110, -80, -100, -110, -10, -90, -10, -30, -30, -50, -70, -150,...
    -50, -130, -110, -90, -70, -40, -60, 0, -20, -100, -80, -90, -110, -30, -10,...
    -70, -50, -40, -60, 0, -20];

chanCoords.y(65:128) = five_by_12_y(65:128);

chanCoords.y(72) = chanCoords.y(72) + 200 
chanCoords.y(122) = chanCoords.y(122) + 400 
chanCoords.y(68) = chanCoords.y(68) + 600 
chanCoords.y(126) = chanCoords.y(126) + 800 


offset = 10*4;
channels = [79,113,78,116,77,115,75,117,80,114,76,118];
for ch = channels
    chanCoords.y(ch) = chanCoords.y(ch) + offset 
end


chanCoords.y(65:128) = chanCoords.y(65:128) - 800;


coords = readtable("D:\github\ripple_heterogeneity\data\electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177 (1).csv");


chanCoords.x(129:end) = coords.x_um_;
chanCoords.y(129:end) = coords.y_um_;

chanCoords.x(129:end) = chanCoords.x(129:end) + 1800

chanCoords.y(129:end) = chanCoords.y(129:end) - 1500
