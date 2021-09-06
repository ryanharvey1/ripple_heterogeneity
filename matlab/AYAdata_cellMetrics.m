clc
%% AYA cell info from excel to cell_metrics

%%
clearvars;
basename = basenameFromBasepath(pwd);
session = sessionTemplate(pwd,'showGUI',false); %
session.channelTags.Bad.channels = [];
save([basename '.session.mat'],'session');

spikes = loadSpikes('basepath',pwd,'format','neurosuite',...
                    'getWaveformsFromDat',false);

%% OR data 
basename = bz_BasenameFromBasepath(pwd);

                
                
%% AYA data
basename = bz_BasenameFromBasepath(pwd);
load([basename '.session.mat']);
load([basename '.spikes.cellinfo.mat']);
load('cell_waveforms.mat');

xfile = dir('*.xlsx');
whowhat = readwhowhat(xfile.name);
  
if numel(spikes.ts) ~= numel(whowhat.id)
    warning('unit n mismatch')
end

% max waveform ch 
for i = 1:numel(spikes.ts)
        for i = 1:numel(cell_metrics.maxWaveformCh)
            try
           cell_metrics.maxWaveformCh1(i) = session.extracellular.electrodeGroups.channels{location(i,1)}(location(i,2));% this is base 1!
           cell_metrics.maxWaveformCh(i) = cell_metrics.maxWaveformCh1(i)-1;
           cell_metrics.waveforms{i} = waveforms(i,:);
            catch
           cell_metrics.maxWaveformCh(i) = NaN;
           cell_metrics.maxWaveformCh1(i) = NaN;
           cell_metrics.waveforms{i} = nan(1,32);
            end
        end 
end

% Pyramidal / interneuron
for i = 1:numel(spikes.ts)
    try
    if whowhat.type(i) == 1
       cell_metrics.putativeCellType{i} = 'Pyramidal Cell';
    elseif whowhat.type(i) == 2
       cell_metrics.putativeCellType{i} = 'Narrow interneuron';
    else
       cell_metrics.putativeCellType{i} = 'Unknown'; 
    end
    catch
       cell_metrics.putativeCellType{i} = 'Unknown';  
    end
end

% region
for i = 1:numel(spikes.ts)
    try
   if whowhat.region(i) == 1
      cell_metrics.brainRegion{i} = 'CA1';
   elseif whowhat.region(i) == 2
      cell_metrics.brainRegion{i} = 'CA2';
   elseif whowhat.region(i) == 3
      cell_metrics.brainRegion{i} = 'CA3';
   elseif whowhat.region(i) == 4
      cell_metrics.brainRegion{i} = 'DG'; 
   else
      cell_metrics.brainRegion{i} = 'Unknown';
   end
    catch
      cell_metrics.brainRegion{i} = 'MEC';  
    end
end
save([basename '.cell_metrics.cellinfo.mat'],'cell_metrics');


%% deep / sup 
load('CA1pos.mat');

if exist('CA1pos','var')
    ii = 0; 
    for i = 1:numel(spikes.ts)
        try
           if whowhat.region(i)== 1 && whowhat.type(i)== 1 && whowhat.layer(i)== 1
               ii = ii +1; 
               if CA1pos(ii,3) == 1
                  cell_metrics.deepSuperficial{i} = 'Deep'; 
               elseif CA1pos(ii,3) == 0
                  cell_metrics.deepSuperficial{i} = 'Superficial'; 
               else
                  cell_metrics.deepSuperficial{i} = 'Unknown'; 
               end
           else
                  cell_metrics.deepSuperficial{i} = 'NA';
           end
        catch
                  cell_metrics.deepSuperficial{i} = 'NA';
        end
    end
    if ii ~= size(CA1pos,1)
        warning('CA1pyr n mismatch')
    end
end

save([basename '.cell_metrics.cellinfo.mat'],'cell_metrics');


