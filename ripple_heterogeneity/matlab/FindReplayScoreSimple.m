function [r,p,a,b] = FindReplayScoreSimple(matrix,varargin)

% FindReplayScore finds the best linear fit in a 2D probability matrix
%
%
% USAGE
%
%    [r,p,a,b] = FindReplayScore(matrix,<options>);
%
%    matrix         
%    threshold      
%    
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'threshold'   considered distance from the line (default = 15 bins)
%     'nShuffles'   default = 500
%     'circular'    for circular-linear data (default = 'on')
%    =========================================================================
%
%   OUTPUT
%
%     r                replay score of matrix
%     p                p-value of replay score
%     a                start position bin of the fitted line
%     b                stop position bin of the fitted line
%
% Copyright (C) 2015 by Ralitsa Todorova and Céline Drieu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Defaults values
nShuffles = 500;
p = [];
circular = 'on';
threshold = 15;

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
end

% Parse parameters
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		builtin('error',['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).']);
	end
	switch(varargin{i}),
		case 'threshold',
			threshold = varargin{i+1};
            if ~isdscalar(threshold,'>0'),
				builtin('error','Incorrect value for property ''threshold'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        case 'nShuffles',
			nShuffles = varargin{i+1};
            if ~isdscalar(nShuffles,'>=0'),
				builtin('error','Incorrect value for property ''nShuffles'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        case 'circular',
            circular = varargin{i+1};
            if ~isstring(circular),
                builtin('error','Incorrect value for property ''circular'' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).');
            end
        otherwise,
			builtin('error',['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindReplayScore">FindReplayScore</a>'' for details).']);
	end
end

% Prepare some default values
nBinsY = size(matrix,1); nBinsX = size(matrix,2); 

%% Get the matrix of sums 
% 'sums' is the score matrix for each possible (x,y) point. 
% The sum of the total score of a replay event would be the sum of the scores of all of its points
% where there are (nBinsX) points, each with its preferred X orientation (among nBinsY choices), 
% where all the points should fit on a line.

matrix = reshape(matrix,nBinsY,[]);

[x,y]= meshgrid((-threshold:1:threshold),1:nBinsY);
indices = mod(x+y-1,nBinsY)+1;
matrixID = zeros(1,1,size(matrix,2));matrixID(:) = 1:size(matrix,2);
indX = reshape(repmat(indices,1,size(matrix,2)),[nBinsY threshold*2+1 size(matrix,2)]);
ind = sub2ind(size(matrix),indX,repmat(matrixID,size(indices)));
sums = squeeze(sum(matrix(ind),2));

%% Get indices describing all possible lines of the sums matrix:

% a and b are the start and end angles (from 1 to nBins) of the linear fit
if strcmp(circular,'on'),
    a = reshape(repmat((1:nBinsY),nBinsY*2-1,1),[],1); % the line can start anywhere
    b = repmat(-(nBinsY-1):(nBinsY-1),1,nBinsY)' + a; % the line has to finish not farther than nBinsY-1 below or above a
else
    a = reshape(repmat((1:nBinsY),nBinsY,1),[],1); % the line can start anywhere from 1 to nBinsY
    b = reshape(repmat((1:nBinsY)',1,nBinsY),[],1); % the line can end anywhere from 1 to nBinsY
end
% make matrix of x-elements 'x'.
% the first column of x is a, and its last column is b
x = [a zeros(length(a),nBinsX-2) b];
% the middle columns go progressively from a to b
for i=2:(nBinsX-1),
    x(:,i) = x(:,1)+(i-1)*(x(:,end)-x(:,1))/(nBinsX-1);
end
% x is rounded so that it designates bin indices
x = mod(round(x)-1,nBinsY)+1;
% y-coordinates
y = repmat(1:nBinsX,size(x,1),1);
indices = sub2ind(size(sums),x,y);

scores = nanmean(sums(indices),2);
[r,ind] = max(scores);
a = a(ind); b = b(ind);

%% Shuffle to get a p-value

if nShuffles>0,
    rShuffled = nan(nShuffles,1);
    aShuffled = nan(nShuffles,1);
    bShuffled = nan(nShuffles,1);
    for i=1:nShuffles,
        shift = round(rand(1,nBinsX)*(nBinsY-1)); % shift columns by a random amount
        shuffledSums = CircularShift(sums,shift);
        [rShuffled(i,1),ind] = max(nanmean(shuffledSums(indices),2));
        aShuffled(i,1) = a(ind); bShuffled(i,1) = b(ind);
    end
    p = sum(rShuffled>r)/nShuffles;    
end