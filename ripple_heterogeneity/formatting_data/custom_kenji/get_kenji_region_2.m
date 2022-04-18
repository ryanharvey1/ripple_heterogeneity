function brain_region = get_kenji_region_2(shankID,basename,dirData)

load([dirData '\ElePosition.mat']);

% pull out shank for this session
shank_region = ElePosition(contains(ElePosition(:,2),basename),6:end);
% make lower case
for i = 1:length(shank_region)
    shank_region{i}=lower(shank_region{i});
end

% add ec to numeric values
for i = 1:length(shank_region)
    if isnumeric(shank_region{i})
        shank_region{i} = ['ec',num2str(shank_region{i})];
    end
end

% locate cell's location by shank region
brain_region = shank_region(shankID);
end