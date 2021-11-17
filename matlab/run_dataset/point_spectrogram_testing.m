
basepath = 'Z:\Data\AYAold\AYA9\day20'
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.ripples.events.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
load(fullfile(basepath,[basename,'.session.mat']))
fs = session.extracellular.sr;
channel = ripples.detectorinfo.detectionparms.Channels(1,1);
lfp = getLFP(channel,'basepath',basepath);
% filter within ripple band
filtered = bz_Filter(lfp,'passband',[100 250]);
filtered.unwrapped = unwrap(filtered.phase);

cell_idx = contains(cell_metrics.putativeCellType,'Pyramidal Cell') & ...
    contains(cell_metrics.brainRegion,'CA1');

% labels = {};
% for i = 1:length(spikes.times)
%     if contains(cell_metrics.deepSuperficial{i},'Superficial')
%         labels{i} = ones(length(spikes.times{i}),1);
%     else
%         labels{i} = ones(length(spikes.times{i}),1)+1;
%     end
% end
% for i = 1:length(spikes.times)
%     labels{i} = zeros(length(spikes.times{i}),1)+i;
% end
for i = 1:length(spikes.times)
    spikes.cycles{i} = interp1(filtered.timestamps,...
        filtered.unwrapped, Restrict(spikes.times{i},ripples.timestamps),...
        'linear')./(2*pi);
end
[ccg,t] = CCG(spikes.cycles(cell_idx),[],'Fs',fs,'binSize',.1,'duration',16,'norm','rate');

labels = cell_metrics.deepSuperficial(cell_idx);
% UID = cell_metrics.UID(cell_idx);

figure;
plot(t,mean(mean(ccg(:,contains(labels,'Deep'),contains(labels,'Deep')),3),2),'DisplayName','deep');
hold on
plot(t,mean(mean(ccg(:,contains(labels,'Superficial'),contains(labels,'Superficial')),3),2),'DisplayName','sup');
plot(t,mean(mean(ccg(:,contains(labels,'Deep'),contains(labels,'Superficial')),3),2),'DisplayName','deep,sup');
legend()
grid('on')
xlabel('ripple cycles')
ylabel('rate')

% figure;
% imagesc(squeeze(ccg(:,1,:)))
% ccg(:,:,1)

figure;
plot(t,mean(mean(ccg(:,contains(labels,'Deep'),contains(labels,'Deep')),3),2));
hold on
x = mean(ccg(:,contains(labels,'Deep'),contains(labels,'Deep')),3)';
plot(t,mean(x)+nanstd(x)/sqrt(size(x,1)),t,mean(x)-nanstd(x)/sqrt(size(x,1)),'color','r')



%%
figure;
imagesc(mean(ccg,3))

cell_metrics.deepSuperficial(cell_idx)


[st,sort_idx] = sort(cat(1,spikes.times{cell_idx}));
label = cat(1,labels{cell_idx});
label = label(sort_idx);

samples = Restrict(st,ripples.timestamps);

samples = Restrict(st,ripples.timestamps);
label = Restrict([st,label],ripples.timestamps);

unwrapped = interp1(filtered.timestamps,unwrap(filtered.phase),samples,'nearest');

[ccg,t] = CCG(unwrapped./(2*pi),label(:,2),...
    'Fs',fs,'binSize',.01,'duration',16,'norm','rate');

unique(label(:,2))
contains(cell_metrics.deepSuperficial{i},'Superficial')

cell_metrics.deepSuperficial(cell_idx)

figure;
plot(t,mean(ccg,3))

figure;
plot(t,ccg(:,1,1))
hold on
plot(t,ccg(:,2,2))
plot(t,ccg(:,1,2))

xlabel('ripple cycles')

figure;
plot(t,ccg(:,1,1)/sum(ccg(:,1,1)))
hold on
plot(t,ccg(:,2,2)/sum(ccg(:,2,2)))
plot(t,ccg(:,1,2)/sum(ccg(:,1,2)))

xlabel('ripple cycles')



% [ccg,t,tau,c] = CCG(unwrapped./(2*pi),label(:,2),'binSize',.1,'duration',16,'mode','ccv');


samples = Restrict([filtered.timestamps,unwrapped],ripples.timestamps);

% restrict spikes to ripples
st = sort(cat(1,spikes.times{cell_idx}));
samples = Restrict(st,ripples.timestamps);
% get unwrapped phase per spike
unwrapped = interp1(filtered.timestamps,unwrap(filtered.phase),samples,'nearest');

cycletime = unwrapped./(2*pi);

[count,cycles] = CountSpikesPerCycle(samples,...
    [filtered.timestamps,filtered.phase]);

spikes_temp.times{1} = cycletime;
PSTH = computePSTH(ripples,spikes_temp);


% compute point spectrogram
% [spectrogram,~,f] = MTPointSpectrogram(cycletime,...
%                                         'frequency',fs,...
%                                         'range',[50, 300],...
%                                         'window',.1,...
%                                         'pad',2,...
%                                         'show','on');

X = nanmean(spectrogram');
figure;
plot(f,X)


% phase_stamps = interp1(filtered.timestamps,filtered.phase,samples,'nearest');
phase_stamps = interp1(filtered.timestamps,unwrapped,samples,'nearest');

phase_stamps = phase_stamps - phase_stamps(1);
window = 500;
interval = zeros(1,ceil(phase_stamps(end)/(2*pi)*window)+window);
interval(1+round(phase_stamps/(2*pi)*window)) = 1;
interval = nanconv(interval,gausswin(120)','edge');
xcorr_spikes2 = xcorr(interval,window);
[~,locs2] = findpeaks(xcorr_spikes2(3*window/2:2*window),'SortStr','descend');


%%
for rip_n = 1:length(ripples.timestamps)
    
    samples = Restrict(st,[ripples.peaks(rip_n,1)-.1,ripples.peaks(rip_n,1)+.1]);
    
    unwrapped = interp1(filtered.timestamps,unwrap(filtered.phase),samples,'nearest');
    

end
%%
[ccg3,t3] = CCG(times2,groups,'binSize',0.0628,'duration',20*pi);

%%

[map,stats] = Map(v,z,varargin)
