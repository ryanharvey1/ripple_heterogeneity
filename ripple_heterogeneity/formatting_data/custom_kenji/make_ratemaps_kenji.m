% make_ratemaps_kenji

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

df = df(contains(df.basepath,'Kenji'),:);

fs = 39.0625;


for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    
    % load position
    if ~exist([basepath,filesep,[basename,'.whl']],'file')
        continue
    end
    
    [x,y,session,spikes] = load_basic_data(basepath,basename,fs);
    
    results = table();
    
    for ep = 1:length(session.epochs)
        time_restrict = [session.epochs{ep}.startTime,session.epochs{ep}.stopTime];
        
        % get ratemaps
        try
            ratemap = get_maps(t,x,y,spikes,3,basepath,time_restrict);
        catch
            continue
        end
        
        for i_unit = 1:length(spikes.times)
            idx = spikes.times{i_unit} >= session.epochs{ep}.startTime &...
                spikes.times{i_unit} <= session.epochs{ep}.stopTime;
            
            si(ep,i_unit) = SpatialInformation('ratemap',ratemap.rateMaps{i_unit}{1},...
                'occupancy',ratemap.occupancy{i_unit}{1},...
                'n_spikes',length(spikes.times{i_unit}(idx)));
        end
        
        % locate place fields
        for i_unit = 1:length(ratemap.rateMaps)
            map.z = ratemap.rateMaps{i_unit}{1};
            map.time = ratemap.occupancy{i_unit}{1};
            map.count = ratemap.countMaps{i_unit}{1};
            map.x = ratemap.params.x;
            map.y = ratemap.params.y;
            stats{i_unit} = MapStats(map);
        end
        
        for i_unit = 1:length(spikes.times)
            specificity(ep,i_unit) = stats{1,i_unit}.specificity;
            n_fields(ep,i_unit) = length(stats{1, i_unit}.size);
        end

        % determine place cell by shuffle
%         shuffle_ic(t,x,y,spikes,basepath,time_restrict,si)
    end
end

function [x,y,session,spikes] = load_basic_data(basepath,basename,fs)
positions = load([basepath,filesep,[basename,'.whl']]);
t = (0:length(positions)-1)'/fs;

% load session and spikes
load([basepath,filesep,[basename,'.session.mat']]);
load([basepath,filesep,[basename,'.spikes.cellinfo.mat']]);

if length(t) ~= length(positions)
    error("sample rate is not correct")
end

% non detects are -1, make then nan
positions(positions == -1) = NaN;

% xy for both leds so make single coords
x = median([positions(:,1),positions(:,3)],2);
y = median([positions(:,2),positions(:,4)],2);
end

function ratemap = get_maps(t,x,y,spikes,binsize,basepath,time_restrict)

idx = t >= time_restrict(1) & t <= time_restrict(2);
temp_t = t(idx);
temp_x = x(idx);
temp_y = y(idx);

% bins
nBins = [round((max(temp_x) - min(temp_x)) / binsize),...
    round((max(temp_y) - min(temp_y)) / binsize)];

% make ratemaps
[ratemap] = firingMapAvg({[temp_t,temp_x,temp_y]},...
    spikes,...
    'nBins',nBins,...
    'basepath',basepath,...
    'saveMat',false,...
    'speedThresh',3);
end

function [p_val,z] = shuffle_ic(t,x,y,spikes,basepath,time_restrict,obs)

temp_spikes = spikes;
for i_unit = 1:length(temp_spikes.times)
    temp_spikes.times{i_unit} =...
        temp_spikes.times{i_unit}(temp_spikes.times{i_unit} >= time_restrict(1) &...
        temp_spikes.times{i_unit} <= time_restrict(2));
end

null = zeros(250,length(temp_spikes.times));
parfor n_shuff = 1:250
    
    ratemap_temp = get_maps(t,x,y,...
        shuff_spike_times(temp_spikes,time_restrict),3,...
        basepath,time_restrict);
    
    null(n_shuff,:) = get_si_for_shuff(temp_spikes,ratemap_temp);
   
end

obs(isnan(obs)) = 0;

p_val = (sum(abs(null) >= abs(obs)) + 1) / (size(null,1) + 1); 
z = (abs(obs) - abs(mean(null))) ./ abs(std(null));

end

function si_shuff = get_si_for_shuff(temp_spikes,ratemap_temp)
    for i_unit = 1:length(temp_spikes.times)
        si_shuff(i_unit) =...
            SpatialInformation('ratemap',ratemap_temp.rateMaps{i_unit}{1},...
            'occupancy',ratemap_temp.occupancy{i_unit}{1},...
            'n_spikes',length(temp_spikes.times{i_unit}));
    end
end

function temp_spikes = shuff_spike_times(temp_spikes,time_restrict)
a = -diff(time_restrict)+20;
b = diff(time_restrict)-20;
r = (b-a).*rand(1,1) + a;
temp_spikes.times = cellfun(@(x) x + r,temp_spikes.times,'un',0);
end

function spatial_information = SpatialInformation(varargin)
% Computes the spatial information of a cell (bits/spike)
%
%
% See parameters below.
%
% See Adrien Peyrache 2008 Methods
%
%  PARAMETERS
%
%   occupation_thresh   Bins that had less than this number of seconds of
%                       occupation are not included in the score. (default 0)
%
%   n_thresh            min number of spikes before value is meaningless(default 50)
%
%   ratemap             pre-computed ratemaps (ncells,nbins)
%
%   occupancy           occupancy maps
%
%   n_spikes            number of spikes

% Modified from TSToolbox by Ryan Harvey 2019

p = inputParser;

p.addParameter('occupation_thresh', 0, @isnumeric);
p.addParameter('n_thresh', 50, @isnumeric);
p.addParameter('ratemap', [], @isnumeric)
p.addParameter('occupancy', [], @isnumeric)
p.addParameter('n_spikes', [], @isnumeric)

p.parse(varargin{:});

occupation_thresh = p.Results.occupation_thresh;
ratemap = p.Results.ratemap;
n_thresh = p.Results.n_thresh;
occupancy = p.Results.occupancy;
n_spikes = p.Results.n_spikes;

% remove low occupancy bins
ratemap(occupancy <= occupation_thresh) = NaN;

% linearize
ratemap = ratemap(:);
occupancy = occupancy(:);

% normalize to probability
occupancy = occupancy/nansum(occupancy);
f = nansum(occupancy.*ratemap);
ratemap = ratemap/f;
ix = ratemap~=0;
SB = (occupancy(ix).*ratemap(ix)).*log2(ratemap(ix));
spatial_information = nansum(SB);

% where spiking was too low, get rid of spatial information score
spatial_information(n_spikes<n_thresh) = NaN;

end