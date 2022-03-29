

% basepath = pwd
% load(fullfile(basepath,[basenameFromBasepath(basepath),'.session.mat']))
% 
% deepSuperficialfromRipple = classification_DeepSuperficial(session)


df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']))
    
    % find polarity reversal shanks
    shanks = [];
    for shank_i = 1:length(deepSuperficialfromRipple.ripple_average)
        [r,c] = size(deepSuperficialfromRipple.ripple_average{shank_i});
        shanks = [shanks;repmat(shank_i,c,1)];
    end
    polarity_reversal = false(length(shanks),1);
    for shank_i = unique(shanks)'
        if length(unique(deepSuperficialfromRipple.channelClass(shanks == shank_i))) > 1
            polarity_reversal(shanks == shank_i) = true;
        end
    end
    
    ripple_average = cell2mat(deepSuperficialfromRipple.ripple_average)';
    ripple_average(polarity_reversal,:)
    ripple_average(~polarity_reversal,:)
    
    df_ripple_info = table();
    df_ripple_info.channel = deepSuperficialfromRipple.channel;
    df_ripple_info.shanks = shanks;
    df_ripple_info.channelClass = deepSuperficialfromRipple.channelClass;
    df_ripple_info.channelDistance = deepSuperficialfromRipple.channelDistance;
    df_ripple_info.polarity_reversal = polarity_reversal;
    df_ripple_info.basepath = repmat(basepath,length(polarity_reversal),1);
    
    [coeff,score,latent,tsquared,explained,mu] = pca(ripple_average);
    df_ripple_info.rip_pc_1 = score(:,1);
    df_ripple_info.rip_pc_2 = score(:,2);
    df_ripple_info.rip_pc_3 = score(:,3);
    df_ripple_info.rip_pc_4 = score(:,4);

    df_ripple_info.rip_pc_1_sqrt = sqrt(df_ripple_info.rip_pc_1+3);
    df_ripple_info.rip_pc_2_sqrt = sqrt(df_ripple_info.rip_pc_2+3);
    df_ripple_info.rip_pc_3_sqrt = sqrt(df_ripple_info.rip_pc_3+3);
    df_ripple_info.rip_pc_4_sqrt = sqrt(df_ripple_info.rip_pc_4+3);
    
    df_ripple_info.rip_pc_1_sqrt = sqrt(df_ripple_info.rip_pc_1+3);
    df_ripple_info.rip_pc_2_sqrt = sqrt(df_ripple_info.rip_pc_2+3);
    df_ripple_info.rip_pc_3_sqrt = sqrt(df_ripple_info.rip_pc_3+3);
    df_ripple_info.rip_pc_4_sqrt = sqrt(df_ripple_info.rip_pc_4+3);
    
    df_ripple_info.rip_pc_1_norm = abs(df_ripple_info.rip_pc_1 - mean(df_ripple_info.rip_pc_1))
    df_ripple_info.rip_pc_2_norm = abs(df_ripple_info.rip_pc_2 - mean(df_ripple_info.rip_pc_2))
    df_ripple_info.rip_pc_3_norm = abs(df_ripple_info.rip_pc_3 - mean(df_ripple_info.rip_pc_3))
    df_ripple_info.rip_pc_4_norm = abs(df_ripple_info.rip_pc_4 - mean(df_ripple_info.rip_pc_4))

    figure;
    histogram(df_ripple_info.rip_pc_1)
    figure;
    histogram(df_ripple_info.rip_pc_1_sqrt)
    
    figure
    scatter(score(:,3),df_ripple_info.channelDistance)
    
    mdl = fitlm(score(:,1:4),df_ripple_info.channelDistance)
    
    mdl = fitlm(df_ripple_info,'channelDistance~rip_pc_1+rip_pc_2+rip_pc_3+rip_pc_4')
    mdl = fitlm(df_ripple_info,'channelDistance~rip_pc_1_norm+rip_pc_2_norm+rip_pc_3_norm+rip_pc_4_norm')
    mdl = fitlm(df_ripple_info,'channelDistance~(rip_pc_1+rip_pc_2+rip_pc_3+rip_pc_4)^2')
    mdl = fitlm(df_ripple_info,'channelDistance~(rip_pc_1_sqrt+rip_pc_2_sqrt+rip_pc_3_sqrt+rip_pc_4_sqrt)^2')
    mdl = fitlm(df_ripple_info,'channelDistance~(rip_pc_1_norm+rip_pc_2_norm+rip_pc_3_norm+rip_pc_4_norm)^2')
    mdl = fitlm(df_ripple_info,'channelDistance~rip_pc_1_norm^2+rip_pc_2_norm^2+rip_pc_3_norm^2+rip_pc_4_norm^2','quadratic')
    mdl = fitlm(df_ripple_info,'channelDistance~rip_pc_1_norm + rip_pc_2_norm^2 + rip_pc_1_norm:rip_pc_4_norm')

    figure;plot(mdl.Fitted,mdl.Residuals.Raw,'.k')
    xlabel('Fitted');ylabel('Residuals.Raw')
    
    figure;plot(mdl.Fitted,mdl.Residuals.Standardized,'.k')
    xlabel('Fitted');ylabel('Residuals.Standardized')
    
    figure;plot(df_ripple_info.rip_pc_1_norm,mdl.Residuals.Standardized,'.k')
    xlabel('rip_pc_1_norm');ylabel('Residuals.Standardized')
    figure;plot(df_ripple_info.rip_pc_2_norm,mdl.Residuals.Standardized,'.k')
    xlabel('rip_pc_2_norm');ylabel('Residuals.Standardized')
    figure;plot(df_ripple_info.rip_pc_3_norm,mdl.Residuals.Standardized,'.k')
    xlabel('rip_pc_3_norm');ylabel('Residuals.Standardized')
    figure;plot(df_ripple_info.rip_pc_4_norm,mdl.Residuals.Standardized,'.k')
    xlabel('rip_pc_4_norm');ylabel('Residuals.Standardized')
    
    figure;
    plot(mdl)
    figure
    plotResiduals(mdl,'probability')
end

function get_data(deepSuperficialfromRipple)
    % find polarity reversal shanks
    shanks = [];
    for shank_i = 1:length(deepSuperficialfromRipple.ripple_average)
        [r,c] = size(deepSuperficialfromRipple.ripple_average);
        shanks = [shanks;repmat(shank_i,c,1)];
    end
    polarity_reversal = false(length(shanks),1);
    for shank_i = unique(shanks)'
        if length(unique(deepSuperficialfromRipple.channelClass(shanks == shank_i))) > 1
            polarity_reversal(shanks == shank_i) = true;
        end
    end
    
    ripple_average = cell2mat(deepSuperficialfromRipple.ripple_average)';
    ripple_average(polarity_reversal,:)
    ripple_average(~polarity_reversal,:)
end