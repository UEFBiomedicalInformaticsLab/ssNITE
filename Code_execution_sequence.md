<table>
  <thead>
    <tr>
      <th align="left">Script Name</th>
      <th align="left">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td colspan="2"><b>Data preparation</b></td>
    </tr>
    <tr>
      <td><a href="Scripts/Build_network/STRING_PPI.R"><code>Scripts/Build_network/STRING_PPI.R</code></a></td>
      <td>Builds a protein-protein interaction network from the STRING database and extracts the largest connected component.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Geneset_libraries/Prepare_geneset_from_databases.R"><code>Scripts/Geneset_libraries/Prepare_geneset_from_databases.R</code></a></td>
      <td>Prepares gene-disease associations by compiling data from numerous databases like OpenTargets, IntOGen, PharmGKB, and others.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Geneset_libraries/KEGGpathway2Gene__human.R"><code>Scripts/Geneset_libraries/KEGGpathway2Gene__human.R</code></a></td>
      <td>Retrieves human KEGG pathway genes and maps them to Ensembl IDs.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Omics_data_prepare/TCGABRCA_data_download.R"><code>Scripts/Omics_data_prepare/TCGABRCA_data_download.R</code></a></td>
      <td>Downloads TCGA-BRCA multi-omics data including RNA-Seq, clinical, proteomics, mutation, CNV, and methylation data.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Omics_data_prepare/TCGABRCA_RNAseq_data_process.R"><code>Scripts/Omics_data_prepare/TCGABRCA_RNAseq_data_process.R</code></a></td>
      <td>Processes RNA-Seq count data by filtering low counts, applying Variance Stabilizing Transformation (VST), and removing outliers.</td>
    </tr>
    <tr>
      <td><a href="Scripts/AssociationsExtract_drug_interactors.R"><code>Scripts/AssociationsExtract_drug_interactors.R</code></a></td>
      <td>Extracts known human protein interactors and targets for small molecular drugs from DrugBank.</td>
    </tr>
    <tr>
      <td><a href="Scripts/AssociationsExtract_adverse_drug_drug_interactions.R"><code>Scripts/AssociationsExtract_adverse_drug_drug_interactions.R</code></a></td>
      <td>Analyzes the drug-drug interactions reported in DrugBank to identify combinations causing severe adverse drug reactions.</td>
    </tr>
    <tr>
      <td colspan="2"><b>Background gene Identification</b></td>
    </tr>
    <tr>
      <td><a href="Scripts/SSN_background_genes/Extract_SSN_background_genes.R"><code>Scripts/SSN_background_genes/Extract_SSN_background_genes.R</code></a></td>
      <td>Short-lists a combined background gene set for network construction using cancer drugs, pathways, and oncogenes, and extracts their STRING sub-networks.</td>
    </tr>
    <tr>
      <td><a href="Scripts/SSN_background_genes/ClusterSamples_SSNbkgGene_expr_TGCABRCA.R"><code>Scripts/SSN_background_genes/ClusterSamples_SSNbkgGene_expr_TGCABRCA.R</code></a></td>
      <td>Checks the clustering of samples via PCA using only the expression levels of the extracted SSN background genes.</td>
    </tr>
    <tr>
      <td colspan="2"><b>Raw single-sample network (SSN) inference</b></td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_ssNITEpy_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_ssNITEpy_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares single-sample networks using the Python-based ssNITEpy random-forest method.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_BONOBO_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_BONOBO_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares single-sample networks using the BONOBO algorithm.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_CSN_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_CSN_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares cell-specific networks (CSN) via a MATLAB backend.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_LIONESS_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_LIONESS_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares single-sample networks using the LIONESS algorithm.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_SWEET_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_SWEET_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares sample-specific-weighted correlation networks (SWEET).</td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_BONOBOsparse_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_BONOBOsparse_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares single-sample networks using the BONOBO algorithm in sparse mode.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_CSNsparse_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_CSNsparse_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares cell-specific networks using the CSN algorithm in sparse mode.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Create_SSN/Create_SWEETsparse_SSN_TCGABRCA.R"><code>Scripts/Create_SSN/Create_SWEETsparse_SSN_TCGABRCA.R</code></a></td>
      <td>Prepares sample-specific networks using the SWEET algorithm in sparse mode.</td>
    </tr>
    <tr>
      <td colspan="2"><b>Network edge refinement & sparsification</b></td>
    </tr>
    <tr>
      <td><a href="Scripts/Edge_weight_stats/EdgeWtStat_SSN_TCGABRCA.R"><code>Scripts/Edge_weight_stats/EdgeWtStat_SSN_TCGABRCA.R</code></a></td>
      <td>Calculates and plots the edge weight statistics (distributions and densities) of the generated SSNs.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Prune_SSN/Prune_SSN_TCGABRCA.R"><code>Scripts/Prune_SSN/Prune_SSN_TCGABRCA.R</code></a></td>
      <td>Prunes the full density SSNs by calculating thresholds based on scale-free topology fitting.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Prune_SSN/Clean_sparse_SSN_TCGABRCA.R"><code>Scripts/Prune_SSN/Clean_sparse_SSN_TCGABRCA.R</code></a></td>
      <td>Cleans inherently sparse SSNs generated by sparse methods by removing zero-weight edges.</td>
    </tr>
    <tr>
      <td colspan="2"><b>Network topology analysis</b></td>
    </tr>
    <tr>
      <td><a href="Scripts/Compute_SSN_centrality/Compute_prunedSSN_centrality_TCGABRCA.R"><code>Scripts/Compute_SSN_centrality/Compute_prunedSSN_centrality_TCGABRCA.R</code></a></td>
      <td>Computes pruned SSNs centralities (degree, closeness, betweenness, etc.) and performs PCA clustering.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Compare_SSN/Compute_prunedSSN_similarity_TCGABRCA.R"><code>Scripts/Compare_SSN/Compute_prunedSSN_similarity_TCGABRCA.R</code></a></td>
      <td>Computes similarity metrics (Jaccard nodes/edges, Adjusted Rand Index) between the pruned SSNs.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Corel_SSNtopol_vs_sampleInfo/Corel_ssnTopol_vs_sampleInfo_TCGABRCA.R"><code>Scripts/Corel_SSNtopol_vs_sampleInfo/Corel_ssnTopol_vs_sampleInfo_TCGABRCA.R</code></a></td>
      <td>Checks the correlation between SSN topological features and sample clinical information.</td>
    </tr>
    <tr>
      <td colspan="2"><b>Prediction of clinical traits</b></td>
    </tr>
    <tr>
      <td><a href="Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/prepare_ML_stratifiedSamples_TCGABRCA.R"><code>Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/prepare_ML_stratifiedSamples_TCGABRCA.R</code></a></td>
      <td>Splits the TCGA-BRCA samples to generate stratified train and test cross-validation splits for ML analysis.</td>
    </tr>
    <tr>
      <td><a href="Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML_geneExpr_vs_sampleInfo_TCGABRCA.R"><code>Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML_geneExpr_vs_sampleInfo_TCGABRCA.R</code></a></td>
      <td>Trains ML models to predict TCGA-BRCA sample clinical traits directly from PCA-reduced gene expression data.</td>
    </tr>
    <tr>
      <td><a href="Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML__ssnTopol_vs_sampleInfo_TCGABRCA.R"><code>Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML__ssnTopol_vs_sampleInfo_TCGABRCA.R</code></a></td>
      <td>Trains ML models to predict TCGA-BRCA sample clinical traits using network-level topological properties.</td>
    </tr>
    <tr>
      <td><a href="Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML_ssnCentralityPC_vs_sampleInfo_TCGABRCA.R"><code>Scripts/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML_ssnCentralityPC_vs_sampleInfo_TCGABRCA.R</code></a></td>
      <td>Trains ML models to predict TCGA-BRCA sample clinical traits from principal components of node/edge centralities.</td>
    </tr>
    <tr>
      <td colspan="2"><b>Cross-method benchmarking & summary visualizations</b></td>
    </tr>
    <tr>
      <td><a href="Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_compare_edge_density.R"><code>Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_compare_edge_density.R</code></a></td>
      <td>Calculates and plots comparative edge weight density distributions across all methods for publication summarization.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_compare_pruneSSN_topology.R"><code>Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_compare_pruneSSN_topology.R</code></a></td>
      <td>Compiles the pruned SSN topological parameters from different methods into comparative boxplot figures.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_prunedSSN_similarity.R"><code>Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_prunedSSN_similarity.R</code></a></td>
      <td>Summarizes the pairwise similarities between SSNs and checks their significance against clinical groups.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_summariseCorel_SSNtopol_vs_sampleInfo.R"><code>Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_summariseCorel_SSNtopol_vs_sampleInfo.R</code></a></td>
      <td>Summarizes the correlation between SSN topological properties and sample information into heatmaps and winner plots.</td>
    </tr>
    <tr>
      <td><a href="Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo.R"><code>Scripts/Compare_X_methods/TCGABRCA/TCGABRCA_summariseMLres_SSNtopol_vs_sampleInfo.R</code></a></td>
      <td>Plots summary grids comparing machine learning prediction accuracies across network topologies, centralities, and gene expression.</td>
    </tr>
  </tbody>
</table>