%% test correlation vs rip particip & n PFs

ses = 'A:\ORproject\Wmaze2\OR15\day10';
basename = bz_BasenameFromBasepath(pwd);

load([basename '.cell_metrics.cellinfo.mat'])
load([basename '.PCsNew.mat']);
load([basename '.SWRunitMetrics.mat']);

for i = 1:numel(cell_metrics.UID)
    if strcmp(cell_metrics.putativeCellType{i},'Pyramidal Cell') && ...
       strcmp(cell_metrics.brainRegion{i},'CA1') && cell_metrics.firingRate(i)<3

       temp = [];
       for j = 1:4
           temp = cat(2,temp,stats{i}{j}.peak);
           cell_metrics.nPF(i,1) = sum(temp>0); 
       end
       
    else
        cell_metrics.nPF(i,1) = NaN;
        SWRunitMetrics.pre.particip(i,1) = NaN;
        SWRunitMetrics.task.particip(i,1) = NaN;
        SWRunitMetrics.post.particip(i,1) = NaN;
    end
end

figure;
subplot(1,3,1);
scatter(cell_metrics.nPF,cell_metrics.firingRate,'.');
xlabel('nPF');ylabel('firing rate');xlim([-1 Inf]);ylim([-0.5 Inf]);
[r,p]=corr(cell_metrics.nPF,cell_metrics.firingRate','rows','complete');
title(['r=' num2str(r,2) ' / p=' num2str(p,2)]);
subplot(1,3,2);
scatter(cell_metrics.nPF,SWRunitMetrics.pre.particip,'.');
xlabel('nPF');ylabel('particip');xlim([-1 Inf]);ylim([-0.1 Inf]);
[r p]=corr(cell_metrics.nPF,SWRunitMetrics.pre.particip,'rows','complete');
title(['r=' num2str(r,2) ' / p=' num2str(p,2)]);
subplot(1,3,3);
scatter(cell_metrics.firingRate,SWRunitMetrics.pre.particip,'.');
xlabel('firing rate');ylabel('particip');xlim([-0.5 Inf]);ylim([-0.1 Inf]);
[r p]=corr(cell_metrics.firingRate',SWRunitMetrics.pre.particip,'rows','complete');
title(['r=' num2str(r,2) ' / p=' num2str(p,2)]);

%%
nPF = cat(1,nPF,cell_metrics.nPF);
FR = cat(1,FR,cell_metrics.firingRate');
particip = cat(1,particip,SWRunitMetrics.pre.particip);

%%
figure;
subplot(1,3,1);
scatter(nPF,FR,'.');
xlabel('nPF');ylabel('firing rate');xlim([-1 Inf]);ylim([-0.5 Inf]);
[r,p]=corr(nPF,FR,'rows','complete');
title(['r=' num2str(r,2) ' / p=' num2str(p,2)]);
subplot(1,3,2);
scatter(nPF,particip,'.');
xlabel('nPF');ylabel('particip');xlim([-1 Inf]);ylim([-0.1 Inf]);
[r p]=corr(nPF,particip,'rows','complete');
title(['r=' num2str(r,2) ' / p=' num2str(p,2)]);
subplot(1,3,3);
scatter(FR,particip,'.');
xlabel('firing rate');ylabel('particip');xlim([-0.5 Inf]);ylim([-0.1 Inf]);
[r p]=corr(FR,particip,'rows','complete');
title(['r=' num2str(r,2) ' / p=' num2str(p,2)]);







