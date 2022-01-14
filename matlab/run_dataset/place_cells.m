df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

% # these datasets don't have multiple mazes in the same session (as far as I know)
not_to_use = contains(df.basepath,'GirardeauG') |... 
        contains(df.basepath,'ORproject') |...
        contains(df.basepath,'OMLproject') |...
        contains(df.basepath,'GrosmarkAD');

df = df(~not_to_use,:);

for i = 1:length(df.basepath)
    placeFieldStats = place_cell_run(df.basepath{i});
    save([basename '.placeFields.cellinfo.mat'],'placeFieldStats');

end

function placeFieldStats = place_cell_run(basepath)
disp(basepath)
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.animal.behavior.mat']))
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

pos{1} = [behavior.time',...
    behavior.position.linearized'];

% Calculate firing maps
[firingMaps] = firingMapAvg(pos,spikes,'basepath',basepath);

% detect place fields
for i = 1:length(firingMaps.rateMaps)
    map.z = firingMaps.rateMaps{i}{1};
    map.time = firingMaps.occupancy{i}{1};
    map.count = firingMaps.countMaps{i}{1};
    map.x = firingMaps.params.x;
    map.y = firingMaps.params.y;
    placeFieldStats{i} = MapStats(map);
end
end