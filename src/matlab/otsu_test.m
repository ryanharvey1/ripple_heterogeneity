
%%
weights = readNPY('C:\Users\Cornell\Downloads\patterns.npy');
weights = weights';
otsuWeights = [];
for assembly = 1:size(weights,2)
    isMember = Otsu(abs(weights(:,assembly)))>1;
    if sum(weights(isMember,assembly)<0)==0 % keep assembly only if there are no negative members
        otsuWeights(:,end+1) = weights(:,assembly).*isMember; % add the weights of the members
    end
end
%%