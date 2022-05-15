

optitrack = optitrack2buzcode('session',session,...
    'basepath', basepath,...
    'basename', 'day10',...
    'filenames','Z:\Data\HMC1\day10\HMC1day10.csv',...
    'unit_normalization',1)

load('day10.DigitalIn.events.mat')

% intanDig = intanDigital2buzcode(session);

OptitrackSync = session.inputs.OptitrackSync.channels; % TTL channel recorded by intan

% Defining timestamps via the TTL-pulses from Optitrack recorded with intan
optitrack.timestamps = digitalIn.timestampsOn{OptitrackSync}(1:numel(optitrack.timestamps));


saveStruct(optitrack,'behavior','session',session);

