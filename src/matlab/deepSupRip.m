function [dsRip] = deepSupRip(spkEventTimes,CA1depth,doPlot,units2include)
% [dsRip] = deepSupRip(spkEventTimes,CA1depth)
% CA1 deep/sup ripple content 

% spkEventTimes = output from bz_getRipSpikes
% CA1depth = distance from middle of pyr layer 

if nargin < 4
   units2include = ones(numel(spkEventTimes.EventAbs),1);
end


    for i = 1:numel(spkEventTimes.EventAbs)
        temp = unique(spkEventTimes.EventAbs{i}(2,:));
        for j = 1:length(temp)
            if units2include(temp(j)) == 0
               temp(j) = NaN;
            end
        end
        dsRip.evtDepth{i,1} = CA1depth(~isnan(temp));
        dsRip.meanDepth(i,1) = nanmean(CA1depth(~isnan(temp)));
        dsRip.deep(i,1) = sum(CA1depth(~isnan(temp))<0)/numel(CA1depth(~isnan(temp)));
        dsRip.sup(i,1) = sum(CA1depth(~isnan(temp))>0)/numel(CA1depth(~isnan(temp)));
        clear temp;
    end

    if doPlot
       figure;
       subplot(1,3,1);
       bar(-5:0.5:5,histc(dsRip.meanDepth,-5:0.5:5));
       title('mean ripple depth')
       subplot(1,3,2);
       bar(0:0.05:1,histc(dsRip.deep,0:0.05:1));
       title('frac deep');
       subplot(1,3,3);
       bar(0:0.05:1,histc(dsRip.sup,0:0.05:1));
       title('frac sup');       
    end
    
end

