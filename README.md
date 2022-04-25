


![Workflow-total_v2](https://user-images.githubusercontent.com/63146629/165044038-53569283-8e15-4345-99c7-d7090511faaa.svg)




# Nicotiana_metabolomics

* **S1_openbabel_conversion.ipynb**	Converts Smiles to Inchi
* **S2_reformater.R**	Clean Databases
* **S3_add_openbabel_info.R**	Merge Databases and drop duplicates
* **S4_run_cfmid_mesocenter.sh** & cfmid_commands.txt	Used to run CFM on HPC Cluster
* **S5_Process_mgf.ipynb**	Used to create composite Spectra from CFM Collision Energies
* **S6_matchms_spec2vec.py**	Used for Database Matching of CFM ID on HPC Cluster
* **S7_matchms_scores_analysis.py**	Used for Database Matching of CFM ID on HPC Cluster
* **S8_MatchMS-v1-cosine-msp.ipynb**	Used for Database Matching of Nicotiana DB
* **S9_MatchMS-v1-cosine.ipynb**	Used for Database Matching of Jassbi
* **S10_Batch-QTOF-sens-v3.xml**	Used for Batch Mode processing of Dataset
* **S11_mgf-rem-redundancy-v4.ipynb**	Remove redundant features
* **S12_Sirius-removev2.ipynb**	Remove redundant IDs from Sirius mgf file
* **S13_run_sirius.sh**	Used to run Sirius on HPC Cluster
* **S14_degree-unsaturation-sirius.ipynb**	Restore Feature ID from Sirius ID and Calculate degree of unsaturation, requires molmass package
* **S15_compound-id-sirius.ipynb**	Restore Feature ID from Sirius ID
* **S16_canopus_consensus_ms2lda.ipynb**	 merge the outputs of all the tools into one big table also get consensus substructures( based on  ms2lda motifs) for insilico-tools and propagate canopus within networks
* **S17_MSLDAmerge-motfs.ipynb**	Get Motifcount based on Presence of Feature
* **S18_MSLDAmerge-motfs-sumall.ipynb**	Get Motifcount based on Presence of Feature within all Tissues
* **S19_canopus_consensus.ipynb**	Script to merge the outputs of all the tools into one big table also get consensus substructures for insilico-tools and propagate canopus within networks
* **S20_phylometabo.ipynb**	Calculate Pairwise Distance Matrix based on Data of Motif Figure and Network Figure
* **S21_phylo.Rmd**	Plot Phylogenies from pairwise Distance Matrix based on APE package
* **S22_sum_molformula_areas.ipynb**	Sum Areas of NANNs based on identical Molecular Formula
* **S23_nann_bubbles.ipynb**	Sum Areas based on Carbon Chain of NANNs, split by hydroxylation or not
* **S24_Networkclustermap.ipynb**	Sum all areas of Networks per Samples
* **S25_ASR-single.Rmd**	Ancestral State reconstruction based on MBASR
* **S26_dbsearch.py**	Script to run cosine score based search on big in-silico db
* **S27_run_db.py**	Script to run cosine score based search on big in-silico db
* **S28_group_for_treemap.ipynb** Used to group and sum canopus classes peak areas
