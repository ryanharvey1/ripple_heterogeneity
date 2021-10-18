basepath = 'A:\Data\AYA9\day20'
load('A:\Data\AYA9\day20\day20.ripples.events.mat')
load('A:\Data\AYA9\day20\day20.cell_metrics.cellinfo.mat')
load('A:\Data\AYA9\day20\day20.spikes.cellinfo.mat')
spikeT = spikes;
%% select events
events=ripples.timestamps(ripples.duration>0.08 & ripples.amplitude>500 ,:);
clear eventsS event
countE=0;
for e = 1:length(events)
    event = events(e,:);
    for i = 1:length(spikeT.times) % colect spk in rip
        temp{i} = Restrict(spikeT.times{i},[event(1) event(end)]);
    end
    countD=0; countS=0;
    for i = 1:length(spikeT.times)
        if ~isempty(temp{i}) && strcmp('Pyramidal Cell',cell_metrics.putativeCellType{i}) && ...
                strcmp('Deep',cell_metrics.deepSuperficial{i})
            countD=countD+1;
        elseif ~isempty(temp{i}) && strcmp('Pyramidal Cell',cell_metrics.putativeCellType{i}) && ...
                strcmp('Superficial',cell_metrics.deepSuperficial{i})
            countS=countS+1;
        end
    end
    if countD > 3 && countS > 3
        countE=countE+1;
        eventsS(countE,:) = events(e,:);
    end
end

%% plot raster
for e = 1:numel(eventsS)
    event = eventsS(e,:);
    plotEventRaster(event,basepath,spikeT,201,'deepSup');
%     suptitle(['rip' num2str(e) 'deepSup']);
    pause;
    
end

%% custom ISIs
spkCA1p  = importSpikes('basepath',basepath,'brainRegion','CA1','cellType','Pyramidal Cell' );
ripSpk = getRipSpikes('basepath',basepath,'events',ripples,'spikes',spkCA1p,'saveMat',false);

isiDD=[];   isiSS=[];   isiDS=[];
for e = 1:numel(ripSpk.EventRel)
    for i = 1:length(ripSpk.EventRel{e})-1
        for j = 1:length(ripSpk.EventRel{e})-1
            if j ~= i && ripSpk.EventRel{e}(2,i) ~= ripSpk.EventRel{e}(2,j)
                if strcmp('Deep',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,i)}) && ...
                        strcmp('Deep',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,j)})
                    isiDD=cat(1,isiDD,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
                elseif strcmp('Superficial',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,i)}) && ...
                        strcmp('Superficial',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,j)})
                    isiSS=cat(1,isiSS,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
                elseif strcmp('Deep',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,i)}) && ...
                        strcmp('Superficial',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,j)}) || ...
                        strcmp('Superficial',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,i)}) && ...
                        strcmp('Deep',cell_metrics.deepSuperficial{ripSpk.EventRel{e}(2,j)})
                    isiDS=cat(1,isiDS,abs(ripSpk.EventRel{e}(1,i)-ripSpk.EventRel{e}(1,j)));
                end
            end
        end
    end
end

meanISI = nan(100000,3);
meanISI(1:numel(isiDS),1) = isiDS*1000;
meanISI(1:numel(isiDD),2) = isiDD*1000;
meanISI(1:numel(isiSS),3) = isiSS*1000;
figure;
boxplot(meanISI);
set(gca,'XTickLabel',{'DS','DD','SS'});

figure;
plot(smooth(histc(isiDS*1000,[1:1:100])/sum(histc(isiDS*1000,[1:1:100])),5,'sgolay'),'k');hold on;
plot(smooth(histc(isiDD*1000,[1:1:100])/sum(histc(isiDD*1000,[1:1:100])),5,'sgolay'),'b');hold on;
plot(smooth(histc(isiSS*1000,[1:1:100])/sum(histc(isiSS*1000,[1:1:100])),5,'sgolay'),'r');hold on;
xlim([0 50]); xlabel('ISI (ms)');

%% ISI and CV2 statistics (from Dan's)- NOT USED
spkCA1p  = importSpikes('basepath',basepath,'brainRegion','CA1','cellType','Pyramidal Cell');

countD=0; countS=0;
for i = 1:length(spkCA1p.UID)
    if strcmp('Deep',cell_metrics.deepSuperficial{spkCA1p.UID(i)})
        countD=countD+1;
        spkCA1d.UID(countD) = spkCA1p.UID(i);
        spkCA1d.times{countD} = spkCA1p.times{i};
    elseif strcmp('Superficial',cell_metrics.deepSuperficial{spkCA1p.UID(i)})
        countS=countS+1;
        spkCA1s.UID(countS) = spkCA1p.UID(i);
        spkCA1s.times{countS} = spkCA1p.times{i};
    end
end

rip.ints.good = ripples.timestamps;
ISIall = ISIStats(spkCA1p,'showfig',true,'ints',rip.ints);
ISId = ISIStats(spkCA1d,'showfig',true,'ints',rip.ints);
ISIs = ISIStats(spkCA1s,'showfig',true,'ints',rip.ints);

meanISI = nan(100,3);
meanISI(1:55,1) = ISIall.summstats.good.meanISI;
meanISI(1:41,2) = ISId.summstats.good.meanISI;
meanISI(1:14,3) = ISIs.summstats.good.meanISI;

figure;
boxplot(meanISI);ylim([0 50]);




