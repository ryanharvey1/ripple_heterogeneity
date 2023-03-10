% read antonio's mat file and write csv for use in python 
load('C:\Users\Cornell\Downloads\testpulsesDS.mat')

df = table();
df.peaks_ratio = [gainD;gainS] * 100;
df.group = [repmat({'Deep'},length(gainD),1); repmat({'Superficial'},length(gainS),1)];

writetable(df,'C:\Users\Cornell\Downloads\testpulsesDS.csv')