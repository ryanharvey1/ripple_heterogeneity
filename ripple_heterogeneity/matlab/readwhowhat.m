function whowhat = readwhowhat(filename)

% USAGE:
%     whatsup(filename)
%
% INPUTS:
%     filename: base name of the .clu and .res file to be read
%     If this was done with aya protocol it should contain a excel file
%     called ratname.xlsx where it shows how good, where was
%     located -region and layer if so, what type and subtype of cell for
%     each. It should also contain a legend with the number labeling for
%     the different characteristic included
%
%
% Aza 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AB1 new
% filebase = ['D:\aza\analysis\data' '\AB1' '\11-16' ];
% ratname='\AB1';
% cluster_group=1:8;
%filename = [filebase '\general' ratname '.xlsx'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Import .xls file with classification of each cells
num = xlsread(filename);

%Organize to output
whowhat.id = num(2:end,1);
whowhat.quality = num(2:end,2);
whowhat.region = num(2:end,3);
whowhat.layer = num(2:end,4); %pyr, rad lac
whowhat.type = num(2:end,5); %pyr int
whowhat.where = num(2:end,6); %shank location
if size(num,2)>6
    whowhat.subtype = num(2:end,6);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Counting neurons
ca1 = whowhat.id(find(whowhat.region==1));
ca2 = whowhat.id(find(whowhat.region==2));
ca3 = whowhat.id(find(whowhat.region==3));
ca4 = whowhat.id(find(whowhat.region==4));
others = whowhat.id(find(whowhat.region==5));
null = whowhat.id(find(whowhat.region==0));

%recounting
if (size(whowhat.id,1)~=(length(ca1)+length(ca2)+length(ca3)+length(ca4)+length(others)+length(null)))
    disp(['Missing someone...dimensions dont match']);
else
    disp(['Everyone rolled call for analysis!']);
    disp('......')
    disp('............')
    disp('.....................')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Counting neurons specifically
ca1pyrlayer = ca1(find(whowhat.layer(ca1)==1));
ca1pyrlayerpyrcell = ca1pyrlayer(find(whowhat.type(ca1pyrlayer)==1));
ca1pyrlayerintcell = ca1pyrlayer(find(whowhat.type(ca1pyrlayer)==2));
ca1orienslayer = ca1(find(whowhat.layer(ca1)==2));
ca1radlayer = ca1(find(whowhat.layer(ca1)==3));
ca1laclayer = ca1(find(whowhat.layer(ca1)==4));
ca1nopyrlayer = cat(1,ca1orienslayer,ca1radlayer,ca1laclayer);


ca2pyrlayer = ca2;
ca2pyrlayerpyrcell = ca2pyrlayer(find(whowhat.type(ca2pyrlayer)==1));
ca2pyrlayerintcell = ca2pyrlayer(find(whowhat.type(ca2pyrlayer)==2));

ca3pyrlayer = ca3;
ca3pyrlayerpyrcell = ca3pyrlayer(find(whowhat.type(ca3pyrlayer)==1));
ca3pyrlayerintcell = ca3pyrlayer(find(whowhat.type(ca3pyrlayer)==2));

ca4pyrlayer = ca4;
ca4pyrlayerpyrcell = ca4pyrlayer(find(whowhat.type(ca4pyrlayer)==1));
ca4pyrlayerintcell = ca4pyrlayer(find(whowhat.type(ca4pyrlayer)==2));

pyrothers = others(find(whowhat.type(others)==1));
intothers = others(find(whowhat.type(others)==2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Overall you have:\n')
disp('......')
disp([num2str(length(ca1)) ' cells in CA1']);
disp([num2str(length(ca2)) ' cells in CA2']);
disp([num2str(length(ca3)) ' cells in CA3']);
disp([num2str(length(ca4)) ' cells in DG']);
disp([num2str(length(others)) ' outside region of interest']);
disp([num2str(length(null)) ' eliminated for very few spikes']);
disp('......')
disp('............')
disp('.....................')
if length(ca1)~=0
    disp('Specifically...in CA1');
    disp([num2str(length(ca1pyrlayerpyrcell)) ' EXC cells in CA1']);
    disp([num2str(length(ca1pyrlayerintcell)) ' INH cells in CA1 pyramidal layer']);
    disp([num2str(length(ca1orienslayer)) ' INH cells in CA1 oriens layer']);
    disp([num2str(length(ca1radlayer)) ' INH cells in CA1 radiatum layer']);
    disp([num2str(length(ca1laclayer)) ' INH cells in CA1 lacunosum-moleculare layer']);
    disp('......')
    disp('............')
    disp('.....................')
end
if length(ca2)~=0
    disp('Specifically...in CA2');
    disp([num2str(length(ca2pyrlayerpyrcell)) ' EXC cells in CA2']);
    disp([num2str(length(ca2pyrlayerintcell)) ' INH cells in CA2']);
    disp('......')
    disp('............')
    disp('.....................')
end
if length(ca3)~=0
    disp('Specifically...in CA3');
    disp([num2str(length(ca3pyrlayerpyrcell)) ' EXC cells in CA3']);
    disp([num2str(length(ca3pyrlayerintcell)) ' INH cells in CA3']);
    disp('......')
    disp('............')
    disp('.....................')
end
if length(ca4)~=0
    disp('Specifically...in DG');
    disp([num2str(length(ca4pyrlayerpyrcell)) ' EXC cells in DG']);
    disp([num2str(length(ca4pyrlayerintcell)) ' INH cells in DG']);
end
end

