import numpy as np
import pandas as pd
from ripple_heterogeneity.assembly import assembly_reactivation, find_sig_assembly


def run(basepath, regions="CA1|PFC|EC1|EC2|EC3|EC4|EC5|MEC", putativeCellType="Pyr"):

    ar = assembly_reactivation.AssemblyReact(
        basepath, brainRegion=regions, putativeCellType=putativeCellType
    )
    # use built in function to load needed data
    ar.load_data()
    # retain pre task post task structure
    ar.restrict_epochs_to_pre_task_post()
    # get weights for the task
    ar.get_weights([ar.epochs[1]])
    # locate sig assemblies and sig members
    patterns, is_member_sig, keep_assembly, is_member = find_sig_assembly.main(
        ar.patterns
    )

    assembly_df = pd.DataFrame()
    assembly_df["patterns"] = patterns.ravel()
    assembly_df["is_member_sig"] = is_member_sig.ravel()
    assembly_df["assembly_n"] = (
        (np.ones_like(patterns).T * np.arange(patterns.shape[0])).T.astype(int).ravel()
    )
    assembly_df["UID"] = np.tile(ar.cell_metrics.UID.values, patterns.shape[0])
    assembly_df["putativeCellType"] = np.tile(
        ar.cell_metrics.putativeCellType.values, patterns.shape[0]
    )
    assembly_df["brainRegion"] = np.tile(
        ar.cell_metrics.brainRegion.values, patterns.shape[0]
    )
    assembly_df["deepSuperficial"] = np.tile(
        ar.cell_metrics.deepSuperficial.values, patterns.shape[0]
    )
    assembly_df["deepSuperficialDistance"] = np.tile(
        ar.cell_metrics.deepSuperficialDistance.values, patterns.shape[0]
    )


    # deep_mec = []
    # deep_pfc = []
    # superficial_mec = []
    # superficial_pfc = []

    # for n in assembly_df.assembly_n.unique():
    #     temp_assembly_df = assembly_df[
    #         (assembly_df.assembly_n == n) & (assembly_df.is_member_sig)
    #     ]
    #     deep_mec.append(
    #         any(temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
    #         & any((temp_assembly_df.deepSuperficial == "Deep"))
    #     )
    #     deep_pfc.append(
    #         any(temp_assembly_df.brainRegion.str.contains("PFC"))
    #         & any((temp_assembly_df.deepSuperficial == "Deep"))
    #     )
    #     superficial_mec.append(
    #         any(temp_assembly_df.brainRegion.str.contains("EC1|EC2|EC3|EC4|EC5|MEC"))
    #         & any((temp_assembly_df.deepSuperficial == "Superficial"))
    #     )
    #     superficial_pfc.append(
    #         any(temp_assembly_df.brainRegion.str.contains("PFC"))
    #         & any((temp_assembly_df.deepSuperficial == "Superficial"))
    #     )
        

    # prop_df = pd.DataFrame()
    # prop_df["prop_cross_region"] = [
    #     np.mean(np.array(superficial_pfc) > 0),
    #     np.mean(np.array(superficial_mec) > 0),
    #     np.mean(np.array(deep_pfc) > 0),
    #     np.mean(np.array(deep_mec) > 0),
    # ]
    # prop_df["labels"] = ["Superficial PFC", "Superficial MEC", "Deep PFC", "Deep MEC"]
    # prop_df["n_assemblies"] = len(superficial_pfc)
    # prop_df['basepath'] = basepath

    # results = {"assembly_df": assembly_df, "prop_df": prop_df}
    # return results

    # target_regions = ["PFC","EC1|EC2|EC3|EC4|EC5|MEC"]
    # for ca1_sub in ["Deep", "Superficial"]:
    #     # iterate over target regions
    #     for region in target_regions:
    #         for n in assembly_df.assembly_n.unique():
    #             temp_assembly_df = assembly_df[
    #                 (assembly_df.assembly_n == n) & (assembly_df.is_member_sig)
    #             ]
