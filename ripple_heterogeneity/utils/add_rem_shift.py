import loading
import numpy as np

def add_rem_shift(df):
    """
    Takes df and adds rem shifting columns
    input: 
        df: data frame with columns [basepath,UID,deepSuperficial]
    output: 
        df with added columns [rem_shift,non_rem_shift,m_rem,m_wake,p_rem,circ_dist,layer_rem_shift]
    """
    df.loc[:,'rem_shift'] = False
    df.loc[:,'non_rem_shift'] = False
    df.loc[:,'m_wake'] = np.nan
    df.loc[:,'m_rem'] = np.nan
    df.loc[:,'p_rem'] = np.nan
    df.loc[:,'p_wake'] = np.nan
    df.loc[:,'circ_dist'] = np.nan

    for basepath in df.basepath.unique():
        df_rem_shift, _ = loading.load_theta_rem_shift(basepath)
        if df_rem_shift.shape[0] == 0:
            continue

        # restrict rem shift df to valid UIDs from df
        df_rem_shift = df_rem_shift[np.in1d(df_rem_shift.UID,df[df.basepath == basepath].UID)]
        
        # get index of current basepath
        idx = df.basepath==basepath
        for df_rem_shift_temp in df_rem_shift.itertuples():
            df.loc[idx & (df.UID == df_rem_shift_temp.UID),"rem_shift"] = df_rem_shift_temp.rem_shift==1
            df.loc[idx & (df.UID == df_rem_shift_temp.UID),"non_rem_shift"] = df_rem_shift_temp.non_rem_shift==1
            df.loc[idx & (df.UID == df_rem_shift_temp.UID),"m_rem"] = df_rem_shift_temp.m_rem
            df.loc[idx & (df.UID == df_rem_shift_temp.UID),"m_wake"] = df_rem_shift_temp.m_wake
            df.loc[idx & (df.UID == df_rem_shift_temp.UID),"p_rem"] = df_rem_shift_temp.p_rem
            df.loc[idx & (df.UID == df_rem_shift_temp.UID),"circ_dist"] = df_rem_shift_temp.circ_dist

    df.loc[df.rem_shift == True,"rem_shift"] = "rem_shifting"
    df.loc[df.rem_shift == False,"rem_shift"] = "unknown"

    df.loc[df.non_rem_shift == True,"non_rem_shift"] = "non_rem_shifting"
    df.loc[df.non_rem_shift == False,"non_rem_shift"] = "unknown" 

    # make layer by rem shift column
    layer_rem_shift = []
    for temp_df in df.itertuples():
        if (temp_df.rem_shift == "rem_shifting") & (temp_df.deepSuperficial == "Deep"):
            layer_rem_shift.append("deep_rem_shift")
        elif (temp_df.non_rem_shift == "non_rem_shifting") & (temp_df.deepSuperficial == "Superficial"):
            layer_rem_shift.append("sup_non_rem_shift")
        else:
            layer_rem_shift.append("unknown")
    df['layer_rem_shift'] = layer_rem_shift
    
    return df
    