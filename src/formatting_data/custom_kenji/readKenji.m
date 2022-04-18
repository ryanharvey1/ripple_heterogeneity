%% to read Kenji Mizuseki's dataset

dirData = 'A:\Data\Kenji\';
dirCWT = 'A:\waveletsDG\';
Fs = 1250; SR = 20000;
% sessions with DG/CA3 and EC: 
sessions ={'ec013.895_902','ec013.906_918','ec013.921_927','ec013.931_942','ec013.944_958','ec013.961_974',...
            'ec013.976_985','ec016.659_674','ec016.682_688','ec016.694_711','ec016.715_735','ec016.740_764',...
            'ec016.769_789','ec016.791_810','ec016.813_831','ec016.835_850','ec016.853_867','ec016.871_889','ec016.893_911',...            
            'ec016.914_932','ec016.934_946','ec016.950_965','ec016.969_986','ec016.1002_1023','ec016.1025_1048',...
            '2006-4-18','2006-4-10','2006-6-12','2006-6-13'};
            % ec016.871_889
            % clu1 corrupted in '2006-6-7'
layers = {'CA3','DG','EC2'}; 

load([dirData '\KenjiData2.mat']);
load([dirData '\ElePosition.mat']);

%%
for ses = 11:18%length(sessions)
    dirfile = [dirData sessions{ses} '\'];
    for i = 1:size(ElePosition,1)
        if strcmp(ElePosition{i,2},sessions{ses})
            ses_code = ElePosition{i,5};
        end
    end
    
    STATES{1} = load([dirfile sessions{ses} '.sts.RUN']);
    STATES{2} = load([dirfile sessions{ses} '.sts.REM']);
    
    [spiket, spikeind,numclus,iEleClu] = ReadCluRes([dirfile sessions{ses}],1:length(dir([dirfile sessions{ses} '.res*'])));
    
    spktSTATES=cell(1,2);spkindSTATES=cell(1,2); 
    for s=1:2
        for i = 1:length(STATES{s}) 
            temp1=spiket((spiket>=STATES{s}(i,1) & spiket<=STATES{s}(i,2)));
            spktSTATES{s} = [spktSTATES{s};temp1];   
            temp2 = spikeind((spiket>=STATES{s}(i,1) & spiket<=STATES{s}(i,2)));
            spkindSTATES{s} = [spkindSTATES{s};temp2];  clear temp1 temp2;
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
    if count ~= numclus
       warning(['inconsistent #cells in ' sessions{ses}]);
    end
    clear count i
    
    ca3p = cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==2,2); % num clus 
    ca3i = cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==2,2);          
    DGe =  cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==3,2);
    DGi =  cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==3,2); 
    ec2p = cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==4,2);
    ec2i = cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==4,2); 
    ec3p = cellSesMap(cleanSes==1 & iCellSes{1}==1 & regionSes==5,2);
    ec3i = cellSesMap(cleanSes==1 & iCellSes{2}==1 & regionSes==5,2);         
   

end



