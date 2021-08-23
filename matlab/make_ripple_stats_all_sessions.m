

% make_ripple_stats_all_sessions
animal = {'Kenji','AB1','AB3','AB4','AYA4','AYA6','AYA7','AYA9','AYA10',...
    'OML5','OML3','OML7','OML8','OML10','OML18','OML19',...
    'Wmaze2\OR15','Wmaze2\OR18','Wmaze3\OR22','Wmaze3\OR21','Wmaze3\OR23',...
    'GrosmarkAD\Cicero','GrosmarkAD\Buddy','GrosmarkAD\Achilles','GrosmarkAD\Gatsby'};

dataDir1 = 'A:\Data\';
dataDir2 = 'A:\OptoMECLEC\';
dataDir3 = 'A:\ORproject\';

if isempty(gcp('nocreate'))
    parpool(10)
end
for a = 1:length(animal)
    disp(animal{a})
    if strncmp('OML',animal{a},3)
        base_path = dataDir2;
    elseif strncmp('Wmaze',animal{a},5)
        base_path = dataDir3;
    else
        base_path = dataDir1;
    end
    files = dir([base_path,...
        animal{a},...
        filesep,'**',filesep,...
        filesep,'**',filesep,...
        '*.ripples.events.mat']);
    
    for f = 1:length(files)
        basepath = files(f).folder;
        basename = bz_BasenameFromBasepath(basepath);
        disp(basepath)
        if exist(fullfile(basepath,[basename,'.ripples.events.mat']),'file')
        
            % load detected ripples 
            load(fullfile(basepath,[basename,'.ripples.events.mat']))

            if isfield(ripples,'duration')
                continue
            end
            
            [ripples,re_run] = add_features(ripples,basepath,basename);
            
            if re_run
                disp([basepath,' something odd with lfp channel in ripples.events'])
                load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
                run_ripple_pipe_kenji(basepath,basename,spikes)
                continue
            end
            save(fullfile(basepath,[basename,'.ripples.events.mat']),'ripples')
        end
                
    end
end
function [ripples,re_run] = add_features(ripples,basepath,basename)
try
    detector = ripples.detectorName;
catch
    detector = ripples.detectorinfo.detectorname;
end
if contains(detector,'bz_DetectSWR')
                
%    Channels = ripples.detectorinfo.detectionparms.Channels(1);
   load(fullfile(basepath,[basename,'.session.mat']))

   samples = bz_GetLFP(ripples.detectorinfo.detectionparms.Channels(1)-1,...
       'basepath',basepath,'basename',basename,'noPrompts',true);

%                lfp_ = LoadBinary(fullfile(basepath,[basename,'.lfp']),...
%                    'frequency',1250,'channels',Channels,...
%                    'nChannels',session.extracellular.nChannels);
%     fs = ripples.detectorinfo.detectionparms.SR;
    passband = ripples.detectorinfo.detectionparms.ripBP;
else
    % extract lfp data from ripples struct
    fs = ripples.detectorinfo.detectionparms.frequency;
    passband = ripples.detectorinfo.detectionparms.passband;

    if isfield(ripples.detectorinfo.detectionparms,'lfp')
        samples.samplingRate = fs;
        samples.data = ripples.detectorinfo.detectionparms.lfp;
        samples.timestamps = [0:1/fs:(length(samples.data)/fs)-1/fs]';
    else
        samples = bz_GetLFP(ripples.detectorinfo.detectionparms.channel(1)-1,...
            'basepath',basepath,'basename',basename,'noPrompts',true);
    end
end
     
% filter signal according to detection parameters
filtered = bz_Filter(samples,'filter','butter','passband',passband,'order', 3);

% Compute instantaneous frequency
unwrapped = unwrap(filtered.phase);
frequency = bz_Diff(medfilt1(unwrapped,12),filtered.timestamps,'smooth',0);
frequency = frequency/(2*pi);

ripples.amplitude = interp1(filtered.timestamps,filtered.amp,ripples.peaks,'linear');
ripples.frequency = interp1(filtered.timestamps,frequency,ripples.peaks,'nearest');
ripples.duration = ripples.timestamps(:,2) - ripples.timestamps(:,1);

% something very odd in some sessions...
%   where peakNormedPower and amplitude are not correlated
%   I'm assuming something bad happened, but can't track down the issue
%   If found, re-run.
if isempty(ripples.peakNormedPower)
    re_run = false;
    return
end
[R,P] = corrcoef(ripples.peakNormedPower,ripples.amplitude);

if length(P)==1
    re_run = true;
    return
end
if P(1,2)>0.001
    re_run = true;
else
    re_run = false;
end
% for event = 1:size(ripples.timestamps,1)
%     idx = samples.timestamps>= ripples.timestamps(event,1) &...
%         samples.timestamps<= ripples.timestamps(event,2);
%     [max_amp(event,1),max_idx] = max(filtered.amp(idx));
%     freq = frequency(idx);
%     max_freq(event,1) = freq(max_idx);
% end

% figure;
% % [2500:3000]*1250;
% interval = [[2500:3000]*1250];
% plot(samples.timestamps(interval),samples.data(interval))
% hold on
% plot(samples.timestamps(interval),filtered.data(interval))
% plot(samples.timestamps(interval),filtered.amp(interval))
% % plot(ripples.peaks,interp1(samples.timestamps,filtered.amp,ripples.peaks),'*r')
% % plot(ripples.peaks,interp1(samples.timestamps,filtered.data,ripples.peaks),'*r')
% % plot(ripples.peaks,interp1(samples.timestamps,double(ripples.detectorinfo.detectionparms.lfp),ripples.peaks),'*r')
% 
% plot(ripples.peaks,ripples.amplitude,'*r')
% plot(ripples.peaks,ripples.peakNormedPower,'*g')
% 
% xlim([min(samples.timestamps(interval)),max(samples.timestamps(interval))])


% rip_n = 2;
% range_ = .2;
% interval = samples.timestamps >= ripples.peaks(rip_n) - range_ &...
%     samples.timestamps <= ripples.peaks(rip_n) + range_;
% figure
% plot(samples.timestamps(interval),samples.data(interval))
% hold on
% plot(samples.timestamps(interval),filtered.data(interval))
% axis tight

end
% basepath = ripples.detectorinfo.detectionparms.basepath;
% basename = bz_BasenameFromBasepath(basepath);
% load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
% 
% run_ripple_pipe_kenji(basepath,basename,spikes)



% load('A:\Data\Kenji\ec013.540_561\ec013.540_561.ripples.events.mat')
% 
% lfp = ripples.detectorinfo.detectionparms.lfp;
% fs = ripples.detectorinfo.detectionparms.frequency;
% timestamps = 0:1/fs:(length(lfp)/fs)-1/fs;
% passband = [100,250];
% filtered = bz_Filter(double(lfp),'filter','butter','passband',passband,'order', 3);
% 
% % figure
% % plot(lfp(1:1250*10))
% % hold on
% % plot(filtered(1:1250*10))
% 
% [maps,data,stats] = bz_RippleStats(filtered,timestamps',ripples,'frequency',fs);
% save()
