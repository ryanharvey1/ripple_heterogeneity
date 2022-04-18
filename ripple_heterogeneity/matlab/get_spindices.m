function spikes = get_spindices(spikes)

spikes.numcells = length(spikes.UID);
for cc = 1:spikes.numcells
    groups{cc}=spikes.UID(cc).*ones(size(spikes.times{cc}));
end
if spikes.numcells>0
    alltimes = cat(1,spikes.times{:}); groups = cat(1,groups{:}); %from cell to array
    [alltimes,sortidx] = sort(alltimes); groups = groups(sortidx); %sort both
    spikes.spindices = [alltimes groups];
end

end