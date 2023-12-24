#!/usr/bin/env python
# coding: utf-8

# # MAP Analysis with MitoCheck Single Cells 
# 
# 
# In this notebook, our goal is to apply the Mean Average Precision (MAP) metric, developed in the [copairs](https://github.com/cytomining/copairs) analysis package.
# We apply this metric to the MitoCheck single-cell dataset, to see the effects of genetic pertubations based on their phenotype. 
# 
# 
# Some links to look at:
# - copairs [repo](https://github.com/cytomining/copairs)
# - MitoCheck github [repo](https://github.com/WayScience/mitocheck_data)
# - MitoCheck zenodo [repo](https://zenodo.org/records/7967386)

# In[2]:


import sys
import logging
import pprint
import pathlib

from copairs.map import run_pipeline
import pandas as pd
import numpy as np
import pprint

# imports src
sys.path.append("../")
from src import utils

# setting up logger
# setting up logger
logging.basicConfig(filename="map_analysis.log",
                    level=logging.DEBUG,
                    format='%(levelname)s:%(asctime)s:%(name)s:%(message)s')


# ## Loading Downloaded Data
# 
# In this section, we load the MitoCheck single-cell datasets, including the training, positive controls, and negative controls. 
# For detailed information about the dataset, please refer to the MitoCheck report mentioned above.
# 
# After downloading the data, we perform formatting by dividing it into two sections. 
# The first section comprises the metadata of each individual cell, while the second section presents all quantified features in a numpy array format.
# 
# This formatting is designed to easily integrate with the copairs `run_pipeline()` function, allowing for easy execution of the analysis.

# In[3]:


# parameters
training_singlecell_data = pathlib.Path("../data/raw/training_data.csv.gz").resolve(strict=True)
pos_control_data = pathlib.Path("../data/raw/normalized_data/positive_control_data.csv.gz").resolve(strict=True)
neg_control_data = pathlib.Path("../data/raw/normalized_data/negative_control_data.csv.gz").resolve(strict=True)


# In[4]:


# # loading in the data into dataframe (~10min loading)
# training_sc_data = pd.read_csv(training_singlecell_data).drop("Unnamed: 0", axis=1)
# pos_control_sc_data = pd.read_csv(pos_control_data)
# neg_control_sc_data = pd.read_csv(neg_control_data)

# # adding the Mitocheck_Phenotypic_Class into the controls  and labels
# pos_control_sc_data.insert(0, "Mitocheck_Phenotypic_Class", "pos_control")
# neg_control_sc_data.insert(0, "Mitocheck_Phenotypic_Class", "neg_control")

# # droping column from trainign data since it does not exist in the controls
# training_sc_data = training_sc_data.drop("Metadata_Object_Outline", axis=1)


# In[5]:


# TODO: delete this later for
training_sc_data = pd.read_parquet("../data/processed/training_sc_data.parquet")
pos_control_sc_data = pd.read_parquet("../data/processed/pos_control_sc_data.parquet").sample(frac=0.01, random_state=42)
neg_control_sc_data = pd.read_parquet("../data/processed/neg_control_sc_data.parquet").sample(frac=0.01, random_state=42)

# adding the Mitocheck_Phenotypic_Class into the controls  and labels
neg_control_sc_data.insert(0, "Mitocheck_Phenotypic_Class", "neg_control")

# adding control labels into the dataset
training_sc_data.insert(1, "Metadata_is_control", 0)
neg_control_sc_data.insert(1, "Metadata_is_control", 1)

# droping column from trainign data since it does not exist in the controls
training_sc_data = training_sc_data.drop("Metadata_Object_Outline", axis=1)


# In[5]:


# parameters for pipeline
random_state = 42
pos_sameby = ["Mitocheck_Phenotypic_Class",]
pos_diffby = ["Metadata_Plate", "Metadata_Well"]

neg_sameby = ["Metadata_Plate"]
neg_diffby = ["Mitocheck_Phenotypic_Class"]
null_size = 100
batch_size = 100


# storing all map results based on postiive and negative controls and feature types
map_results_neg_cp = []
map_results_neg_dp = []
map_results_neg_cp_dp = []

# running process
# for loop selects one single phenotype
# then splits the data into metadata and raw feature values
# two different groups that contains 3 splits caused by the types of features
# applie the copairs pipeline
for phenotype in list(training_sc_data["Mitocheck_Phenotypic_Class"].unique()):

    # select training dataset based on phenotype
    selected_training = training_sc_data.loc[training_sc_data["Mitocheck_Phenotypic_Class"] == phenotype]

    # concatenate to positive and negative control
    training_w_neg = pd.concat([selected_training, neg_control_sc_data])

    # spliting metadata and raw feature values
    logging.info("splitting data set into metadata and raw feature values")
    negative_training_cp_meta, negative_training_cp_feats = utils.split_data(training_w_neg, dataset="CP")
    negative_training_dp_meta, negative_training_dp_feats = utils.split_data(training_w_neg, dataset="DP")
    negative_training_cp_dp_meta, negative_training_cp_dp_feats = utils.split_data(training_w_neg, dataset="CP_and_DP")


    # execute pipeline on negative control with trianing dataset with cp features
    logging.info(f"Running pipeline on CP features using {phenotype} phenotype")
    cp_negative_training_result = run_pipeline(meta=negative_training_cp_meta,
                                            feats=negative_training_cp_feats,
                                            pos_sameby=pos_sameby,
                                            pos_diffby=pos_diffby,
                                            neg_sameby=neg_sameby,
                                            neg_diffby=neg_diffby,
                                            batch_size=batch_size,
                                            null_size=null_size)

    # adding shuffle label column
    cp_negative_training_result["shuffled"] = False

    # writing out into csv file
    cp_negative_training_result.to_csv(f"{phenotype}_cp_MAP.csv", index=False)

    # append to list
    map_results_neg_cp.append(cp_negative_training_result)

    # execute pipeline on negative control with trianing dataset with dp features
    logging.info(f"Running pipeline on DP features using {phenotype} phenotype")
    dp_negative_training_result = run_pipeline(meta=negative_training_dp_meta,
                                            feats=negative_training_dp_feats,
                                            pos_sameby=pos_sameby,
                                            pos_diffby=pos_diffby,
                                            neg_sameby=neg_sameby,
                                            neg_diffby=neg_diffby,
                                            batch_size=batch_size,
                                            null_size=null_size)

    # adding shuffle label column
    dp_negative_training_result["shuffled"] = False

    # writing out into csv file
    dp_negative_training_result.to_csv(f"{phenotype}_dp_MAP.csv", index=False)

    map_results_neg_dp.append(dp_negative_training_result)

    # execute pipeline on negative control with trianing dataset with cp_dp features
    logging.info(f"Running pipeline on CP and DP features using {phenotype} phenotype")
    cp_dp_negative_training_result = run_pipeline(meta=negative_training_cp_dp_meta,
                                            feats=negative_training_cp_dp_feats,
                                            pos_sameby=pos_sameby,
                                            pos_diffby=pos_diffby,
                                            neg_sameby=neg_sameby,
                                            neg_diffby=neg_diffby,
                                            batch_size=batch_size,
                                            null_size=null_size)

    # adding shuffle label column
    cp_dp_negative_training_result["shuffled"] = False

    # writing out into csv file
    cp_dp_negative_training_result.to_csv(f"{phenotype}_cp_dp_MAP.csv", index=False)

    map_results_neg_cp_dp.append(cp_dp_negative_training_result)


    logging.info("Running MAP Pipeline with shuffled data")
    # Below, we are running the same test, but we are shuffling the phenotypes
    logging.info("Shuffling data based on the Mitocheck_Phenotypic_Class (phenotype) labels")
    shuffled_labels = np.random.permutation(training_w_neg["Mitocheck_Phenotypic_Class"])
    training_w_neg["Mitocheck_Phenotypic_Class"] = shuffled_labels

    # split the shuffled dataset
    # spliting metadata and raw feature values
    logging.info("splitting shuffled data set into metadata and raw feature values")
    shuffled_negative_training_cp_meta, shuffled_negative_training_cp_feats = utils.split_data(training_w_neg, dataset="CP")
    shuffled_negative_training_dp_meta, shuffled_negative_training_dp_feats = utils.split_data(training_w_neg, dataset="DP")
    shuffled_negative_training_cp_dp_meta, shuffled_negative_training_cp_dp_feats = utils.split_data(training_w_neg, dataset="CP_and_DP")


    # execute pipeline on negative control with trianing dataset with cp features
    logging.info(f"Running pipeline on CP features using {phenotype} phenotype, data is shuffled")
    shuffled_cp_negative_training_result = run_pipeline(meta=shuffled_negative_training_cp_meta,
                                            feats=shuffled_negative_training_cp_feats,
                                            pos_sameby=pos_sameby,
                                            pos_diffby=pos_diffby,
                                            neg_sameby=neg_sameby,
                                            neg_diffby=neg_diffby,
                                            batch_size=batch_size,
                                            null_size=null_size)

    # adding shuffle label column
    shuffled_cp_negative_training_result["shuffled"] = True

    # writing out into csv file
    shuffled_cp_negative_training_result.to_csv(f"{phenotype}_cp_shuffled_MAP.csv", index=False)

    # append to list
    map_results_neg_cp.append(shuffled_cp_negative_training_result)

    # execute pipeline on negative control with trianing dataset with dp features
    logging.info(f"Running pipeline on DP features using {phenotype} phenotype, data is shuffled")
    shuffled_dp_negative_training_result = run_pipeline(meta=shuffled_negative_training_dp_meta,
                                            feats=shuffled_negative_training_dp_feats,
                                            pos_sameby=pos_sameby,
                                            pos_diffby=pos_diffby,
                                            neg_sameby=neg_sameby,
                                            neg_diffby=neg_diffby,
                                            batch_size=batch_size,
                                            null_size=null_size)

    # adding shuffle label column
    shuffled_dp_negative_training_result["shuffled"] = True

    # writing out into csv file
    shuffled_dp_negative_training_result.to_csv(f"{phenotype}_dp_shuffled_MAP.csv", index=False)

    map_results_neg_dp.append(shuffled_dp_negative_training_result)

    # execute pipeline on negative control with trianing dataset with cp_dp features
    logging.info(f"Running pipeline on CP and DP features using {phenotype} phenotype, data is shuffled")
    shuffled_cp_dp_negative_training_result = run_pipeline(meta=shuffled_negative_training_cp_dp_meta,
                                            feats=shuffled_negative_training_cp_dp_feats,
                                            pos_sameby=pos_sameby,
                                            pos_diffby=pos_diffby,
                                            neg_sameby=neg_sameby,
                                            neg_diffby=neg_diffby,
                                            batch_size=batch_size,
                                            null_size=null_size)

    # adding shuffle label column
    shuffled_cp_dp_negative_training_result["shuffled"] = True

    # writing out into csv file
    shuffled_cp_dp_negative_training_result.to_csv(f"{phenotype}_cp_dp_shuffled_MAP.csv", index=False)

    map_results_neg_cp_dp.append(shuffled_cp_dp_negative_training_result)

