
% add_familiarity_to_session
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'Kenji'),:);

py_file = "D:\github\ripple_heterogeneity\ripple_heterogeneity\utils\compress_repeated_epochs.py"
disp(pyversion)
pathToCode = fileparts(which(py_file));
if count(py.sys.path,pathToCode)==0
    % adds the code path to the python path
    insert(py.sys.path,int32(0),pathToCode) 
end

for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    
    disp(basepath)
    
    load(fullfile(basepath,[basename,'.session.mat']))
    
    epoch_df = [];
    for ep = 1:length(session.epochs)
        epoch_df{ep} = struct2cell(session.epochs{ep})';
    end
    epoch_df = cell2table(vertcat(epoch_df{:}),"VariableNames",fields(session.epochs{1}));
    epoch_df = table2struct(epoch_df,"ToScalar",true)
    test = py.dict(epoch_df)
    
    

    pyrunfile("py_file epoch_df","epoch_df")
    pyOut = py.compress_repeated_epochs.compress_repeated_epochs(vertcat(epoch_df{:}));
    %     save(fullfile(basepath,[basename,'.session.mat']),'session')
end

% function epoch_df = compress_repeated_epochs(epoch_df)
% %     """
% %     Compress back to back sleep epoch
% %     """
%     match = zeros([epoch_df.environment.shape[0]])
%     match[match == 0] = np.nan
%     for i,ep in enumerate(epoch_df.environment[:-1]):
%         if np.isnan(match[i])
% %             # find match in current and next epoch
%             if (ep == epoch_df.environment.iloc[i+1]) & (ep == 'sleep'):
%                 match[i:i+2] = i
% %                 # given match, see if there are more matches
%                 for match_i in np.arange(1,epoch_df.environment[:-1].shape[0])
%                     if i+1+match_i == epoch_df.environment.shape[0]
%                         break
%                     end
%                     if ep == epoch_df.environment.iloc[i+1+match_i]
%
%                         match[i:i+1+match_i+1] = i
%                     else
%                         break
%                     end
%                 end
%             end
%         end
%     for i in range(len(match))
%         if np.isnan(match[i])
%             match[i] = (i+1)*2000 %# make nans large numbers that are unlikely to be real epoch
%         end
%     end
% %     # iter through each epoch indicator to get start and stop
%     results = pd.DataFrame()
%     no_nan_match = match[~np.isnan(match)]
%     for m in pd.unique(no_nan_match):
%         temp_dict = {'name': epoch_df[match==m].name.iloc[0],
%                     'startTime':epoch_df[match==m].startTime.iloc[0],
%                     'stopTime':epoch_df[match==m].stopTime.iloc[-1],
%                     'environment':epoch_df[match==m].environment.iloc[0],
%                     'behavioralParadigm':epoch_df[match==m].behavioralParadigm.iloc[0],
%                     }
%         temp_df = pd.DataFrame.from_dict(temp_dict,orient='index').T
%
%         results = results.append(temp_df,ignore_index=True)
%     end
%     end
% end