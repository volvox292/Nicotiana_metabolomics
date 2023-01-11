# Evolutionary metabolomics of specialized metabolism diversification in the genus *Nicotiana* highlights allopolyploidy-mediated innovations in N-acylnornicotine metabolism

preprint available at: https://doi.org/10.1101/2022.09.12.507566

supplementatry data at : [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6536010.svg)](https://doi.org/10.5281/zenodo.6536010)

```bash {bash, echo=T, eval=F}
# Running QCxmS on the singularity container (ubuntu.sif) with the qcxms.py script 
# —> singularity exec ubuntu.sif python qcxms.py -o $outdir -nt $nr_trajectories -i $inchi  

# Sending it 
module load singularity
srun -p public singularity exec ubuntu.sif python qcxms.py -o mol1 -nt 5 -i "InChI=1S/C23H38N2O3/c1-2-3-5-11-20(26)12-6-4-7-13-21(27)17-23(28)25-16-9-14-22(25)19-10-8-15-24-18-19/h8,10,15,18,20-22,26-27H,2-7,9,11-14,16-17H2,1H3"
```


## Workflow:




![Workflow-total_v2](https://user-images.githubusercontent.com/63146629/165044114-0ca1595a-4b74-4bb4-9e8b-f83161431060.png)







## Scripts used in this publication:

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
* **S29_alpha_diversity.ipynb** Calculate alpha diversity based on shannon entropy
* **S30_Vegan_calculations.Rmd** NMDS using vegan package
* **S31_cosine_distance_sp_canopus.ipynb** Calculate distances between species and CANOPUS classes
