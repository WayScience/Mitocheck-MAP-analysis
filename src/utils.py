"""
Contains utitlity functions to import into the anatysis


`split_data()` was develoepd by @roshankern:
https://github.com/WayScience/mitocheck_data/blob/63f37859d993b8de25fefe1cb8a3aac421c3e08a/utils/load_utils.py#L84
"""
import pandas as pd

def split_data(pycytominer_output: pd.DataFrame, dataset: str = "CP_and_DP"):
    """
    split pycytominer output to metadata dataframe and np array of feature values

    Parameters
    ----------
    pycytominer_output : pd.DataFrame
        dataframe with pycytominer output
    dataset : str, optional
        which dataset features to split,
        can be "CP" or "DP" or by default "CP_and_DP"

    Returns
    -------
    pd.Dataframe, np.ndarray
        metadata dataframe, feature values
    
    Credit: 
    """
    all_cols = pycytominer_output.columns.tolist()

    # get DP,CP, or both features from all columns depending on desired dataset
    if dataset == "CP":
        feature_cols = [col for col in all_cols if "CP__" in col]
    elif dataset == "DP":
        feature_cols = [col for col in all_cols if "DP__" in col]
    elif dataset == "CP_and_DP":
        feature_cols = [col for col in all_cols if "P__" in col]

    # metadata columns is all columns except feature columns
    metadata_cols = [col for col in all_cols if "P__" not in col]

    metadata_dataframe = pycytominer_output[metadata_cols]
    feature_data = pycytominer_output[feature_cols].values

    return metadata_dataframe, feature_data