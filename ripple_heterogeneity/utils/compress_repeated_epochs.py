import numpy as np
import pandas as pd


def main(epoch_df:pd.DataFrame, epoch_name:str="sleep") -> pd.DataFrame:
    """
    Compress back to back epochs with name epoch_name
    Input:
        epoch_df: pandas dataframe with epochs
        epoch_name: name of epochs to compress
    Output:
        epoch_df: pandas dataframe with epochs
    """
    match = np.zeros([epoch_df.environment.shape[0]])
    match[match == 0] = np.nan
    for i, ep in enumerate(epoch_df.environment[:-1]):
        if np.isnan(match[i]):
            # find match in current and next epoch
            if (ep == epoch_df.environment.iloc[i + 1]) & (ep == epoch_name):
                match[i : i + 2] = i
                # given match, see if there are more matches
                for match_i in np.arange(1, epoch_df.environment[:-1].shape[0]):
                    if i + 1 + match_i == epoch_df.environment.shape[0]:
                        break
                    if ep == epoch_df.environment.iloc[i + 1 + match_i]:

                        match[i : i + 1 + match_i + 1] = i
                    else:
                        break

    for i in range(len(match)):
        if np.isnan(match[i]):
            match[i] = (
                i + 1
            ) * 2000  # make nans large numbers that are unlikely to be real epoch

    # iter through each epoch indicator to get start and stop
    results = pd.DataFrame()
    no_nan_match = match[~np.isnan(match)]
    for m in pd.unique(no_nan_match):
        temp_dict = {}
        for item in epoch_df.keys():
            temp_dict[item] = epoch_df[match == m][item].iloc[0]

        temp_dict["startTime"] = epoch_df[match == m].startTime.min()
        temp_dict["stopTime"] = epoch_df[match == m].stopTime.max()

        temp_df = pd.DataFrame.from_dict(temp_dict, orient="index").T

        results = pd.concat([results, temp_df], ignore_index=True)
    return results
