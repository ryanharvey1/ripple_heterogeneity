% pulses_to_manipulation

% Data container for manipulation data. A MATLAB struct
%   manipulationName stored in a .mat file: basename.eventName.manipulation.mat with the following fields:
% 
% timestamps: Px2 matrix with intervals for the P events in seconds.
% peaks: Event time for the peak of each events in seconds (Px1).
% amplitude: amplitude of each event (Px1).
% amplitudeUnits: specify the units of the amplitude vector.
% eventID: numeric ID for classifying various event types (Px1).
% eventIDlabels: cell array with labels for classifying various event types defined in stimID (cell array, Px1).
% eventIDbinary: boolean specifying if eventID should be read as binary values (default: false).
% center: center time-point of event (in seconds; calculated from timestamps; Px1).
% duration: duration of event (in seconds; calculated from timestamps; Px1).
% detectorinfo: info about how the events were detected.
df = table();
df.basepath = {'Z:\Data\Can\OML22\day6',...
    'Z:\Data\Can\OML22\day7',...
    'Z:\Data\Can\OML22\day8',...
    'Z:\Data\Can\OML22\day19',...
    'Z:\Data\Can\OML22\day20'}';

for i = 1:length(df.basepath)
   run(df.basepath{i}) 
end

function run(basepath)

basename = basenameFromBasepath(basepath);

load(fullfile(basepath,[basename,'.pulses.events.mat']))
load(fullfile(basepath,[basename,'.session.mat']))

if ~isfield(pulses,'timestamps')
    pulses.timestamps = pulses.intsPeriods;
end
if ~isfield(pulses,'amplitude')
    pulses.amplitude = nan(length(pulses.timestamps),1);
end
if ~isfield(pulses,'eventGroupID')
    pulses.eventGroupID = ones(length(pulses.timestamps),1);
end
if ~isfield(pulses,'duration')
    pulses.duration = pulses.timestamps(:,2) - pulses.timestamps(:,1);
end
optoStim.timestamps = pulses.timestamps;
optoStim.peaks = median(pulses.timestamps,2);
optoStim.amplitude = pulses.amplitude;
optoStim.amplitudeUnits = 'pulse_respect_baseline';
optoStim.eventID = pulses.eventGroupID;
optoStim.eventIDlabels = {'laser_on'};
optoStim.eventIDbinary = true;
optoStim.center = median(pulses.timestamps,2);
optoStim.duration = pulses.duration;
optoStim.detectorinfo = 'getAnalogPulses_pulses_to_manipulation';

saveStruct(optoStim,'manipulation','session',session);
end