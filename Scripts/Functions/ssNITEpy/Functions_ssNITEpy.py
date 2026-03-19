# Functions regarding ssNITEpy


# Define the function to build patient specific network

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from collections import defaultdict
from operator import itemgetter
from multiprocessing import Pool
from copy import copy
import numpy as np
import pandas as pd
from scipy import sparse
import warnings
import os

def ssNITEpy(
    expr_data, 
    gene_names = None, 
    target_genes = 'all', 
    regulators = 'all', 
    rf_args = {'n_estimators':100, 'random_state':0}, 
    n_trees = 50, 
    max_links = 'all',
    n_folds = 5, 
    n_repeats = 3, 
    frequency_threshold = 0.5, 
    n_workers = 1, 
    return_adjacency = False, 
    decision_path_importance = False, 
    link_save_file = None
):
    """
    Builds patient-specific regulatory networks using K-Fold cross-validation with repeated folds,
    prioritizing robust links based on their frequency across folds and repetitions.
    
    Parameters:
    - expr_data (pd.DataFrame): Gene expression data (samples x genes).
    - gene_names (list): List of gene names corresponding to the columns of expr_data.
    - target_genes (list or 'all'): List of targetd genes or 'all' to use all genes.
    - regulators (list or 'all'): List of candidate regulator genes or 'all' to use all genes.
    - rf_args (dict): Arguments passed to the RandomForestRegressor for each gene.
    - n_trees (int): Number of top trees to consider for each patient.
    - max_links (int or 'all'): Number of top links to retain in the output. Use 'all' for all links.
    - n_folds (int): Number of folds for cross-validation.
    - n_repeats (int): Number of repetitions of K-Fold cross-validation.
    - frequency_threshold (float): Minimum fraction of folds/repetitions in which a link must appear.
    - n_workers (int): Number of parallel workers to use. 
    - return_adjacency (bool) : whether to return adjacency matrix instead of link list pd.DataFrame
    - decision_path_importance (bool) : whether to use patient decision paths for importance instead of tree subset
    - link_save_file (str) : if defined and return_adjacency is False, writes each target's link lists into files ending in link_save_file
    Returns:
    - dict with keys "patient_networks" and "error_summary"
    """
    if gene_names is None:
        gene_names = expr_data.columns.to_numpy()
    if not isinstance(gene_names, np.ndarray):
        gene_names = np.array(gene_names)
    
    if not isinstance(target_genes, (list, np.ndarray)) and target_genes != 'all':
        raise ValueError("'target_genes' must be a list of gene names or 'all'.")
    else:
        target_genes = gene_names
    
    
    link_save_file = link_save_file + "ssNITEpy_edges_tmp/"
    print(link_save_file)
    if not os.path.exists(link_save_file): 
        os.makedirs(link_save_file, exist_ok = True)

    
    if n_workers > 1:
        patient_networks = []
        patient_errors = []

        with Pool(processes = n_workers, maxtasksperchild = 1) as mp:
            for gene_patient_networks, gene_patient_errors in mp.starmap(
                build_patient_specific_networks_with_folds, 
                [
                    (expr_data, 
                     gene_names, 
                     [i], 
                     regulators, 
                     rf_args, 
                     n_trees, 
                     n_folds, 
                     n_repeats, 
                     max_links,
                     frequency_threshold, 
                     decision_path_importance, 
                     return_adjacency, 
                     link_save_file) 
                    for i in target_genes
                ]
            ):
                patient_networks.append(gene_patient_networks)
                patient_errors.append(gene_patient_errors)
        if not link_save_file:
            if return_adjacency:
                patient_networks = np.concatenate(patient_networks, axis = 2)
            else:
                patient_networks = pd.concat(patient_networks, axis = 0)
        patient_errors = np.concatenate(patient_errors, axis = 1)
    else:
        patient_networks, patient_errors = build_patient_specific_networks_with_folds(
            expr_data = expr_data, 
            gene_names = gene_names, 
            target_genes = target_genes, 
            regulators = regulators, 
            rf_args = rf_args, 
            n_trees = n_trees, 
            n_folds = n_folds, 
            n_repeats = n_repeats, 
            max_links = max_links,
            frequency_threshold = frequency_threshold, 
            decision_path_importance = decision_path_importance, 
            return_adjacency = return_adjacency, 
            link_save_file = link_save_file
        )
    
    patient_errors = summarise_patient_errors(patient_errors) 
    patient_errors = patient_errors.set_index(expr_data.index)
    
    # return patient_networks, patient_errors
    return {"patient_networks" : patient_networks, "error_summary" : patient_errors}


def build_patient_specific_networks_with_folds(
    expr_data, 
    gene_names=None, 
    target_genes='all',
    regulators='all', 
    rf_args = {'n_estimators':100, 'random_state':0}, 
    n_trees=50, 
    n_folds=5, 
    n_repeats=3, 
    max_links='all',
    frequency_threshold=0.5, 
    decision_path_importance = False, 
    return_adjacency = False, 
    link_save_file = None
):
    if not isinstance(expr_data, pd.DataFrame):
        raise ValueError("expr_data must be a pandas DataFrame.")
    
    if gene_names is None:
        gene_names = expr_data.columns.to_numpy()
    elif not isinstance(gene_names, np.ndarray):
        gene_names = np.array(gene_names)
    
    if isinstance(target_genes, (list, np.ndarray)):
        pass
    elif target_genes == 'all':
        target_genes = gene_names
    else:
        raise ValueError("'target_genes' must be a list of gene names or 'all'.")
    
    if isinstance(regulators, (list, np.ndarray)):
        n_regs = len(regulators)
    elif regulators == 'all':
        n_regs = len(gene_names)
    else:
        raise ValueError("'regulators' must be a list of gene names or 'all'.")
    
    # Initialize storage for errors per patient
    patient_errors = np.zeros((expr_data.shape[0], len(target_genes), n_repeats), dtype='float32')
    
    patient_networks = []
    # Iterate over each target gene
    for target_gene_idx, target_gene in enumerate(target_genes):
        
        # Initialize patient-specific link storage for target gene
        target_patient_links = np.zeros((expr_data.shape[0], n_regs, n_repeats), dtype='float32')
        
        # Repeated K-Fold Cross-Validation
        for repeat in range(n_repeats):
            kf = KFold(n_splits=n_folds, shuffle=True, random_state=repeat)
            
            for fold, (train_idx, test_idx) in enumerate(kf.split(expr_data)):
                # Split the data
                train_data = expr_data.iloc[train_idx]
                test_data = expr_data.iloc[test_idx]
                
                # Target variable: expression levels of the target gene
                y_train = train_data[target_gene]
                y_test = test_data[target_gene]
                
                # Candidate regulators
                if regulators == 'all':
                    candidate_regulators_mask = gene_names != target_gene
                else:
                    candidate_regulators = [gene for gene in regulators if gene != target_gene]
                    candidate_regulators_mask = [gene in candidate_regulators for gene in gene_names]
                    candidate_regulators_mask = np.array(candidate_regulators_mask)
                
                # Features: expression levels of the candidate regulators
                X_train = train_data.loc[:,candidate_regulators_mask]
                X_test = test_data.loc[:,candidate_regulators_mask]
                
                # Train a Random Forest
                rf = RandomForestRegressor(**rf_args)
                rf.fit(X_train, y_train)
                
                if decision_path_importance:
                    node_indicator, _ = rf.decision_path(X_test)
                    tree_node_features = np.concatenate([i.tree_.feature for i in rf.estimators_], axis = 0)
                    for patient_test_idx in np.arange(X_test.shape[0], dtype = 'int'):
                        _, patient_node_ptr = node_indicator[patient_test_idx].nonzero()
                        patient_features = tree_node_features[patient_node_ptr]
                        patient_features = patient_features[patient_features >= 0] # don't count leaf nodes
                        gene_idx, gene_count = np.unique(patient_features, return_counts = True)
                        masked_gene_idx = np.argwhere(candidate_regulators_mask)[gene_idx, 0]
                        gene_importance = gene_count / np.sum(gene_count)
                        patient_idx = test_idx[patient_test_idx]
                        # Insert importances to target_patient_links
                        target_patient_links[patient_idx, masked_gene_idx, repeat] = gene_importance
                    y_test_pred = rf.predict(X_test)
                    target_error = np.abs(y_test_pred - y_test.to_numpy())
                    patient_errors[test_idx, target_gene_idx, repeat] = target_error
                else:
                    # Extract the trained trees
                    trees = rf.estimators_
                    importance_matrix = np.array([tree.feature_importances_ for tree in trees])
                    
                    # Get top trees for each sample in the test set
                    errors_matrix = np.zeros((X_test.shape[0], len(trees)))
                    for tree_idx, tree in enumerate(trees):
                        predictions = tree.predict(X_test.to_numpy())
                        errors = np.abs(y_test.to_numpy() - predictions)
                        errors_matrix[:, tree_idx] = errors
                    
                    top_trees_matrix = np.argsort(errors_matrix, axis = 1)[:,:n_trees]
                    patient_errors[test_idx, target_gene_idx, repeat] = np.mean(np.take_along_axis(errors_matrix, top_trees_matrix, axis = 1), axis = 1)
                    
                    patient_importances = np.take_along_axis(
                        np.expand_dims(importance_matrix, axis = 0), 
                        np.expand_dims(top_trees_matrix, axis = 2), 
                        axis = 1
                    )
                    patient_importances_mean = np.mean(patient_importances, axis = 1)
                    # Find index in full array and insert
                    insert_idx = np.ix_(
                        test_idx, 
                        candidate_regulators_mask, 
                        [repeat]
                    )
                    target_patient_links[insert_idx] = patient_importances_mean[:, :, np.newaxis]
        
        # De-duplicate here to reduce memory foot-print
        target_patient_networks = deduplicate_and_prioritize_links(
            target_patient_links = target_patient_links, 
            frequency_threshold = frequency_threshold, 
            max_links = max_links
        )
        
        if not return_adjacency:
            # convert to DataFrame of edge lists
            nnan_ind = np.logical_not(np.isnan(target_patient_networks))
            link_idx = np.argwhere(nnan_ind)
            link_val = target_patient_networks[nnan_ind]
            id_vec = expr_data.index[link_idx[:,0]]
            reg_gene_vec = gene_names[link_idx[:,1]]
            #target_gene_vec = gene_names[link_idx[:,2]]
            target_patient_networks = pd.DataFrame(
                {
                    'ID' : id_vec, 
                    'regulator_gene' : reg_gene_vec, 
                    'target_gene' : target_gene, #target_gene_vec, 
                    'importance' : link_val
                }
            )
        if link_save_file and not return_adjacency:
            target_patient_networks.to_parquet(
                f"{link_save_file}ssNITEpy_edges_{target_gene}.parquet.lz4", 
                compression = "lz4",
                index = False
            )
            patient_networks.append(f"{link_save_file}ssNITEpy_edges_{target_gene}.parquet.lz4")
        else:
            patient_networks.append(target_patient_networks)
    
    if return_adjacency:
        patient_networks = np.transpose(np.array(patient_networks), axes = (1,2,0))
    elif not link_save_file:
        patient_networks = pd.concat(patient_networks, axis = 0)
    
    return patient_networks, patient_errors

def deduplicate_and_prioritize_links(
    target_patient_links, 
    frequency_threshold, 
    max_links
):
    # Deduplicate and prioritize robust links
    zero_mask = target_patient_links == 0
    target_patient_links[zero_mask] = np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        patient_networks = np.nanmean(target_patient_links, axis = 2)
    if frequency_threshold > 0.:
        freq = np.mean(zero_mask, axis = 2)
        freq_mask = freq > (1. - frequency_threshold)
        patient_networks[freq_mask] = np.nan
    if isinstance(max_links, int) and max_links > 0:
        bottom_regs = np.argsort(-patient_networks, axis = 1)
        for patient_idx in range(patient_networks.shape[0]):
            for target_idx in range(patient_networks.shape[2]):
                bottom_regs_ij = bottom_regs[patient_idx,max_links:,target_idx]
                patient_networks[patient_idx,bottom_regs_ij,target_idx] = np.nan
    
    return patient_networks

def summarise_patient_errors(
    patient_errors
):
    # Summarize error distribution for each patient
    error_summary_df = pd.DataFrame({
        'mean_error' : np.mean(patient_errors, axis = (1,2)),
        'median_error' : np.median(patient_errors, axis = (1,2)),
        'std_error' : np.std(patient_errors, axis = (1,2)),
        'min_error' : np.min(patient_errors, axis = (1,2)),
        'max_error' : np.max(patient_errors, axis = (1,2))
    })
    return error_summary_df



#####



# Define the function to plot the error density 

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde

def plot_cluster_density(data, column, min_cluster_size=4, filename=None):
    """
    Plot density plots for each cluster based on a specified column.

    Parameters:
    - data (pd.DataFrame): DataFrame containing error statistics and cluster labels.
    - column (str): Column name to plot density for.
    - min_cluster_size (int): Minimum number of samples required in a cluster to compute density.
    """
    plt.figure(figsize=(10, 6))
    x = np.linspace(data[column].min() - 0.02, data[column].max() + 0.02, 100)

    for cluster in sorted(data["Cluster"].unique()):
        cluster_data = data[data["Cluster"] == cluster][column]

        # Skip clusters with fewer samples than the threshold
        if len(cluster_data) < min_cluster_size:
            print(f"Skipping Cluster {cluster + 1} due to insufficient data (only {len(cluster_data)} samples)")
            continue

        # Compute density using Gaussian KDE
        kde = gaussian_kde(cluster_data)
        y = kde(x)

        plt.plot(x, y, label=f"Cluster {cluster + 1}")

    plt.title(f"Density Plot for {column.capitalize()} by Cluster")
    plt.xlabel(column.capitalize())
    plt.ylabel("Density")
    plt.legend(title="Clusters", loc="upper right")
    plt.grid(True)
    plt.tight_layout()
    if filename is not None:
      plt.savefig(filename)
    else:
        plt.show()
    plt.close() 
