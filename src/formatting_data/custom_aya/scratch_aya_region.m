% scratch_aya_region
df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');
basepath = unique(df(contains(df.animal,'AYA'),:).basepath);

for i = 1:length(basepath)
    disp(basepath{i})
    add_region(basepath{i})
end


cell_metrics = loadCellMetricsBatch('basepaths',basepath');
cell_metrics = CellExplorer('metrics',cell_metrics);

function add_region(basepath)
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
    if isfield(cell_metrics,'CA1depth') &&...
            isfield(cell_metrics,'brainRegion') &&...
            isfield(cell_metrics,'brainRegion')
        return
    end
    try
        cell_metrics2 = load(fullfile(basepath,[basename,'.cell_metrics.cellinfo_old.mat']));
    catch
        try
            cell_metrics2 = load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.legacy.mat']));
        catch
            cell_metrics2 = load(fullfile(basepath,[basename,'.cell_metrics.cellinfo1.mat']));
        end
    end
    cell_metrics.brainRegion = cell_metrics2.cell_metrics.brainRegion;
    cell_metrics.CA1depth(1,:) = cell_metrics2.cell_metrics.CA1depth;
    
    save(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']),'cell_metrics')
end