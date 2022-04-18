

function [group,threshold,em] = Otsu(vector)
 
% The Otsu method for splitting data into two groups.
% This is somewhat equivalent to kmeans(vector,2), but while the kmeans implementation
% finds a local minimum and may therefore produce different results each time,
% the Otsu implementation is guaranteed to find the best division every time.
%
% EXAMPLE USAGE:
% [~,~,weights] = ActivityTemplatesICA(spikes,'bins',intervals);
% otsuWeights = [];
% for assembly = 1:size(weights,2)
%     isMember = Otsu(abs(weights(:,assembly)))>1;
%     if sum(weights(isMember,assembly)<0)==0 % keep assembly only if there are no negative members
%         otsuWeights(:,end+1) = weights(:,assembly).*isMember; % add the weights of the members
%     end
% end
 
sorted = sort(vector);
intraClassVariance = nan(size(vector));
n = length(vector);
for i=1:n-1
    p = (i)/n; p0 = 1-p;
    intraClassVariance(i) = p*var(sorted(1:i),1)+ p0*var(sorted(i+1:end),1);
end
[minIntraVariance,idx] = min(intraClassVariance);
threshold = sorted(idx);
group = (vector > threshold)+1;
 
em = 1 - (minIntraVariance/var(vector,1)); % em = effectiveness metric
 
%% Alternative code (from wikipedia)
% [ histogramCounts,ht] = hist(vector,length(vector));
% total = sum(histogramCounts); % total number of pixels in the image
% top = total;
% sumB = 0;
% wB = 0;
% maximum = 0.0;
% sum1 = dot(0:top-1, histogramCounts);
% for ii = 1:top
%     wF = total - wB;
%     if wB > 0 && wF > 0
%         mF = (sum1 - sumB) / wF;
%         val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);
%         if ( val >= maximum )
% 	  level = ii;
% 	  maximum = val;
%         end
%     end
%     wB = wB + histogramCounts(ii);
%     sumB = sumB + (ii-1) * histogramCounts(ii);
% end
% 
% threshold = ht(level);
% group = (vector > threshold)+1;
end
