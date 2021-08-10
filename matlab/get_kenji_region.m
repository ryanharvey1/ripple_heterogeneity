function brain_region=get_kenji_region(varargin)
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'basename',[],@isstr); 
addParameter(p,'dirData','A:\Data\Kenji\',@isstr); 
addParameter(p,'check_cell_count',true,@logical); 

% Parsing inputs
parse(p,varargin{:})
basepath = p.Results.basepath;
basename = p.Results.basename;
dirData = p.Results.dirData;
check_cell_count = p.Results.check_cell_count;

if isempty(basename)
    parts = strsplit(basepath,filesep);
    basename = parts{end};
end

dirfile = [dirData, basename, filesep];

load([dirData '\KenjiData2.mat']);
load([dirData '\ElePosition.mat']);

for i = 1:size(ElePosition,1)
    if strcmp(ElePosition{i,2},basename)
        ses_code = ElePosition{i,5};
    end
end

% Get cell types
cellSesMap=[];cleanSes=[];regionSes=[];iCellSes=cell(1,2); count=0;
for i = 1:size(PyrIntMap.Map,1)
    if PyrIntMap.Map(i,18) == ses_code  
       count = count+1;
       cellSesMap(count,:) = PyrIntMap.Map(i,:);
       cleanSes(count,:) = Clean(i);
       iCellSes{1}(count,:) = iCell{1}(i); % 1=pyr
       iCellSes{2}(count,:) = iCell{2}(i);% 1= int
       regionSes(count,:) = Region(i);
    end    
end

if check_cell_count
    [~,~,numclus,~] = ReadCluRes(...
                                [dirfile,basename],...
                                1:length(dir([dirfile,basename,'.res*']))...
                                );
    if count ~= numclus
       error(['inconsistent #cells in ' basename]);
    end
end

regions = {'ca1','ca3','dg','ec2','ec3','ec4','ec5','ecq'};
% brain_region = regions(cellSesMap(:,19));

brain_region = regions(regionSes);

% ca3p = cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==2,2);
% ca3i = cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==2,2);          
% DGe =  cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==3,2);
% DGi =  cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==3,2); 
% ec2p = cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==4,2);
% ec2i = cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==4,2); 
% ec3p = cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==5,2);
% ec3i = cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==5,2);         

% brain_region = repmat({'unknown'},1,length(regionSes));
% brain_region(ca3p) = {'CA3p'};
% brain_region(ca3i) = {'CA3i'};
% brain_region(DGe) = {'DGe'};
% brain_region(DGi) = {'DGi'};
% brain_region(ec2p) = {'EC2p'};
% brain_region(ec2i) = {'EC2i'};
% brain_region(ec3p) = {'EC3p'};
% brain_region(ec3i) = {'EC3i'};





end