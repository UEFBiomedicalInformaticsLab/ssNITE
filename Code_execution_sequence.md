<h2 style="font-family: sans-serif; color: #333; margin-bottom: 10px;">Code execution sequence</h2>

<table style="width: 100%; border-collapse: collapse; font-family: sans-serif;">

<thead>
<tr style="background-color: #f2f2f2;">
<th style="width: 40%; padding: 10px; border: 1px solid #ddd; text-align: left;">Script Name</th>
<th style="width: 60%; padding: 10px; border: 1px solid #ddd; text-align: left;">Description</th>
</tr>
</thead>

<tbody>

<tr style="background-color: #e9ecef;">
<td colspan="2" style="padding: 12px; border: 1px solid #ddd; text-align: center;"><strong>Data preparation</strong></td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Build_network/STRING_PPI.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Builds a protein-protein interaction network from the STRING database and extracts the largest connected component.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Geneset_libraries/Prepare_geneset_from_databases.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Prepares gene-disease associations by compiling data from numerous databases like OpenTargets, IntOGen, PharmGKB, and others.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Geneset_libraries/KEGGpathway2Gene__human.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Retrieves human KEGG pathway genes and maps them to Ensembl IDs.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Omics_data_prepare/TCGABRCA_data_download.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Downloads TCGA-BRCA multi-omics data including RNA-Seq, clinical, proteomics, mutation, CNV, and methylation data.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Omics_data_prepare/TCGABRCA_RNAseq_data_process.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Processes RNA-Seq count data by filtering low counts, applying Variance Stabilizing Transformation (VST), and removing outliers.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/AssociationsExtract_drug_interactors.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Extracts known human protein interactors and targets for small molecular drugs from DrugBank.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/AssociationsExtract_adverse_drug_drug_interactions.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Analyzes the drug-drug interactions reported in DrugBank to identify combinations causing severe adverse drug reactions.</td>
</tr>

<tr style="background-color: #e9ecef;">
<td colspan="2" style="padding: 12px; border: 1px solid #ddd; text-align: center;"><strong>Background gene Identification</strong></td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/SSN_background_genes/Extract_SSN_background_genes.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Short-lists a combined background gene set for network construction using cancer drugs, pathways, and oncogenes.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/SSN_background_genes/ClusterSamples_SSNbkgGene_expr_TGCABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Checks the clustering of samples via PCA using only the expression levels of the extracted SSN background genes.</td>
</tr>

<tr style="background-color: #e9ecef;">
<td colspan="2" style="padding: 12px; border: 1px solid #ddd; text-align: center;"><strong>Raw single-sample network (SSN) inference</strong></td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Create_SSN/Create_ssNITEpy_SSN_TCGABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Prepares single-sample networks using the Python-based ssNITEpy random-forest method.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Create_SSN/Create_BONOBO_SSN_TCGABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Prepares single-sample networks using the BONOBO algorithm.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Create_SSN/Create_CSN_SSN_TCGABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Prepares cell-specific networks (CSN) via a MATLAB backend.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Create_SSN/Create_LIONESS_SSN_TCGABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Prepares single-sample networks using the LIONESS algorithm.</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Create_SSN/Create_SWEET_SSN_TCGABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Prepares sample-specific-weighted correlation networks (SWEET).</td>
</tr>

<tr style="background-color: #e9ecef;">
<td colspan="2" style="padding: 12px; border: 1px solid #ddd; text-align: center;"><strong>Network edge refinement & sparsification</strong></td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Edge_weight_stats/EdgeWtStat_SSN_TCGABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Calculates and plots the edge weight statistics (distributions and densities).</td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Prune_SSN/Prune_SSN_TCGABRCA.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Prunes full density SSNs based on scale-free topology fitting.</td>
</tr>

<tr style="background-color: #e9ecef;">
<td colspan="2" style="padding: 12px; border: 1px solid #ddd; text-align: center;"><strong>Prediction of clinical traits</strong></td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/ML_SSN.../ML__ssnTopol_vs_sampleInfo.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Trains ML models to predict traits using network-level topological properties.</td>
</tr>

<tr style="background-color: #e9ecef;">
<td colspan="2" style="padding: 12px; border: 1px solid #ddd; text-align: center;"><strong>Cross-method benchmarking & summary</strong></td>
</tr>
<tr>
<td style="padding: 8px; border: 1px solid #ddd;"><code>Scripts/Compare_X_methods/.../summariseMLres.R</code></td>
<td style="padding: 8px; border: 1px solid #ddd;">Plots summary grids comparing prediction accuracies across methods.</td>
</tr>

</tbody>
</table>