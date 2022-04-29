df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'OR') | contains(df.basepath,'OML'),:);

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'OR'),:);
general_behavior_file('basepath',df.basepath)


for i = 1:length(df.basepath)
    disp(df.basepath{i})
    extract_tracking(df.basepath{i})
end
function extract_tracking(basepath)
% basepath = 'Z:\Data\OMLproject\OML7\day7'
basename = basenameFromBasepath(basepath);
animal_id = animalFromBasepath(basepath);

trackingFiles = dir([basepath,filesep,animal_id,basename,'*.csv']);

% bad_idx = OR15day11.RigidBody7 > 0.005;
% x = OR15day11.RigidBody4;
% y = OR15day11.RigidBody6;
% x(bad_idx) = NaN;
% y(bad_idx) = NaN;
% 
% figure;plot(x,y)

% Get csv file locations
nFiles = length(trackingFiles);

% Load merge point info to correct times
mergeFile = checkFile('basepath',basepath,'fileType','.MergePoints.events.mat');
load([mergeFile.folder,filesep,mergeFile.name]);

% Load data from each file with proper time
clear trackData;
for fileIdx = 1:nFiles
    tempFile = trackingFiles(fileIdx);
    trackData.file(fileIdx) = readOptitrackCSV(tempFile.name,'basepath',tempFile.folder);
    tempFolder = tempFile(fileIdx).folder;
    subDirIdx = contains(extractBefore(tempFile(fileIdx).name,'.csv'),MergePoints.foldernames);
    t0 = MergePoints.timestamps(subDirIdx,1);
    trackData.file(fileIdx).timestamps = trackData.file(fileIdx).timestamps + t0;
end


% Stick data together in order
data = vertcat(trackData.file(:).data);
timestamps = vertcat(trackData.file(:).timestamps);

tracking = struct();
tracking.timestamps = timestamps;
tracking.frameCount = data(:,1);

tracking.orientation.rx = data(:,3);
tracking.orientation.ry = data(:,4);
tracking.orientation.rz = data(:,5);
tracking.orientation.rw = data(:,6);

tracking.position.x = data(:,7);
tracking.position.y = data(:,8);
tracking.position.z = data(:,9);

if size(data,2)==10
    tracking.errorPerMarker = data(:,10);
end
end