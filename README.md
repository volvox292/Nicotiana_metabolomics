# Nicotiana_metabolomics

* S1	openbabel_conversion.ipynb	Converts Smiles to Inchi
* S2	reformater.R	Clean Databases
* S3	add_openbabel_info.R	Merge Databases and drop duplicates
* S4	run_cfmid_mesocenter.sh & cfmid_commands.txt	Used to run CFM on HPC Cluster
* S5	Process_mgf.ipynb	Used to create composite Spectra from CFM Collision Energies
* S6	matchms_spec2vec.py	Used for Database Matching of CFM ID on HPC Cluster
* S7	matchms_scores_analysis.py	Used for Database Matching of CFM ID on HPC Cluster
* S8	MatchMS-v1-cousine-msp.ipynb	Used for Database Matching of Nicotiana DB
* S9	MatchMS-v1-cousine.ipynb	Used for Database Matching of Jassbi
* S10	Batch-QTOF-sens-v3.xml	Used for Batch Mode processing of Dataset
* S11	mgf-rem-redundancy-v4.ipynb	Remove redundant features
* S12	Sirius-removev2.ipynb	Remove redundant IDs from Sirius mgf file
* S13	run_sirius.sh	Used to run Sirius on HPC Cluster
* S14	degree-unsaturation-sirius.ipynb	Restore Feature ID from Sirius ID and Calculate degree of unsaturation, requires molmass package
* S15	compound-id-sirius.ipynb	Restore Feature ID from Sirius ID
* S16	canopus_consensus_ms2lda.ipynb	 merge the outputs of all the tools into one big table also get consensus substructures( based on  ms2lda motifs) for insilico-tools and propagate canopus within networks
* S17	MSLDAmerge-motfs.ipynb	Get Motifcount based on Presence of Feature
* S18	MSLDAmerge-motfs-sumall.ipynb	Get Motifcount based on Presence of Feature within all Tissues
* S19	canopus_consensus.ipynb	Script to merge the outputs of all the tools into one big table also get consensus substructures for insilico-tools and propagate canopus within networks
* S20	phylometabo.ipynb	Calculate Pairwise Distance Matrix based on Data of Motif Figure and Network Figure
* S21	phylo.Rmd	Plot Phylogenies from pairwise Distance Matrix based on APE package
* S22	sum_molformula_areas.ipynb	Sum Areas of NANNs based on identical Molecular Formula
* S23	nann_bubbles.ipynb	Sum Areas based on Carbon Chain of NANNs, split by hydroxylation or not
* S24	Networkclustermap.ipynb	Sum all areas of Networks per Samples
* S25	ASR-single.Rmd	Ancestral State reconstruction based on MBASR
* S26	dbsearch.py	Script to run cosine score based search on big in-silico db
* S27	run_db.py	Script to run cosine score based search on big in-silico db
