
Z:\Data\FujisawaS\EE\EE0627fm
idx = behavior.timestamps > session.epochs{1, 2}.startTime & behavior.timestamps < session.epochs{1, 2}.stopTime

figure;
plot(behavior.position.x(idx),behavior.position.y(idx))

figure;
plot(behavior.timestamps,behavior.position.x)

intervals = []
y_lim = ylim
for ep = 1:length(session.epochs)
%     intervals(ep,:) = [session.epochs{1, ep}.startTime,session.epochs{1, ep}.stopTime];
    PlotIntervals([session.epochs{1, ep}.startTime,session.epochs{1, ep}.stopTime],'color',rand(1,3))
    text(session.epochs{1, ep}.startTime,y_lim(2),session.epochs{1, ep}.name)
    text(session.epochs{1, ep}.startTime,y_lim(1),session.epochs{1, ep}.environment)
end

for ep = 1:length(session.epochs)
%     intervals(ep,:) = [session.epochs{1, ep}.startTime,session.epochs{1, ep}.stopTime];
[session.epochs{1, ep}.stopTime- session.epochs{1, ep}.startTime]/60
end
% PlotIntervals(intervals,'alphaValue',0.5)

figure
for ep = 1:length(SessionNP)
%     intervals(ep,:) = [session.epochs{1, ep}.startTime,session.epochs{1, ep}.stopTime];
    PlotIntervals([SessionNP(ep,1),SessionNP(ep,2)],'color',rand(1,3))
end
SessionNP

behavior.timestamps = behavior.timestamps + session.epochs{1, 4}.startTime  


behavior.timestamps(behavior.timestamps > session.epochs{1, 4}.stopTime) =...
    behavior.timestamps(behavior.timestamps > session.epochs{1, 4}.stopTime) + session.epochs{1, 9}.startTime  

behavior.timestamps = behavior.timestamps + session.epochs{1, 4}.startTime  
