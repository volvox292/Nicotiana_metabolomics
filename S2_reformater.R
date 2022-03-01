library(tidyverse)
library(data.table)


# Information -------------------------------------------------------------

# We will load and clean each database here.
# To be inserted in the end to the sqlite database, each database needs an unique ID, SMILES, InChI, InChIKey.
# We will use openbabel to retrieve missing SMILES, InChI and InChIKey.


# 2 columns databases substances ----------------------------------------
# Load all the substances DB formats (14 ATM)
substances <- list.files("Ressources/0_Raw_db_files", "*substances.csv", full = T)
names(substances) <- basename(substances)
substances.dt <- rbindlist(lapply(substances, fread), idcol = "Source")
setnames(substances.dt, c("zinc_id", "smiles"), c("ID", "SMILES"))

# See duplicates
substances.dt[duplicated(substances.dt[, c("ID", "SMILES")]), ]

# remove duplicates 
substances.nodup.dt <- substances.dt[!duplicated(substances.dt[, c("ID", "SMILES")]), ]

# Example of duplicated entry
substances.dt[ID == "ZINC000004098020"]
# Should now have one entry
substances.nodup.dt[ID == "ZINC000004098020"]


fwrite(substances.nodup.dt, "Ressources/1_Clean_db_files/SUBS.csv")

# One column file for openbabel
fwrite(substances.nodup.dt[, c("SMILES")], "Ressources/2_Openbabel_conversion/SUBS.smiles", col.names = F)

# We can convert SMILES to INCHI but its very slow :( with 
#library(ChemmineOB)
#substances.nodup.dt[, InChI := convertFormat("SMI", "INCHI", smiles), by = smiles]
#convertFormatFile(from = "SMI", to = "INCHI", fromFile = "Ressources/Openbabel_conversion/SUBS.smiles", toFile = "Ressources/Openbabel_conversion/SUBS.inchi")

# UNPD --------------------------------------------------------------------
# Has smiles, inchi and inchikey information
unpd.dt <- fread("Ressources/0_Raw_db_files/UNPD.csv")
colnames(unpd.dt)
unpd.dt[, Source := "UNPD"]
setnames(unpd.dt, c("inchik", "PMA_ID"), c("InChIKey", "ID"))
unpd.subset.dt <- unpd.dt[, c("Source", "ID", "SMILES", "InChI", "InChIKey")]

# Remove duplicated InChI
unpd.subset.nodup.dt <- unpd.subset.dt[!duplicated(unpd.subset.dt[, c("InChI")]), ]
fwrite(unpd.subset.nodup.dt, "Ressources/1_Clean_db_files/UNPD.csv")

# Knapsack ----------------------------------------------------------------
# Only has InChI
knapsack <- fread("Ressources/0_Raw_db_files/list_knapsack.txt")
knapsack[, Source := "knapsack"]
colnames(knapsack)
setnames(knapsack, c("C_ID", "INCHI-CODE"), c("ID","InChI"))

knapsack.dt <- knapsack[, c("Source", "ID", "InChI")]
knapsack.nodup.dt <- knapsack.dt[!duplicated(knapsack.dt[, c("InChI")]), ]
fwrite(knapsack.nodup.dt, "Ressources/1_Clean_db_files/KNAP.csv")

fwrite(knapsack.nodup.dt[, c("InChI")], "Ressources/2_Openbabel_conversion/KNAP.inchi", col.names = F)

# Mexican Natural Products ------------------------------------------------
# has only SMILES
mexicanNP <- fread("Ressources/0_Raw_db_files/apps_database_csv_BIOFACQUIM-Mexican-Natural-Products.csv")
mexicanNP[, Source := "Mexican Natural Products"]
mexicanNP.dt <- mexicanNP[, c("Source", "ID", "SMILES")]
mexicanNP.nodup.dt <- mexicanNP.dt[!duplicated(mexicanNP.dt[, c("SMILES")]), ]

fwrite(mexicanNP.nodup.dt, "Ressources/1_Clean_db_files/MENP.csv")

# for openbabel
fwrite(mexicanNP.nodup.dt[, c("SMILES")], "MENP.smiles", col.names = F)


# LOTUS_DB ----------------------------------------------------------------
# Only has SMILES
lotusDB <- fread("Ressources/0_Raw_db_files/LOTUS_DB.smi", header = F, col.names = c("SMILES", "lotus_id"))
lotusDB[, Source := "LotusDB"]
colnames(lotusDB)
setnames(lotusDB, "lotus_id", "ID")
lotusDB.nodup.dt <- lotusDB[!duplicated(lotusDB[, c("SMILES")])]
fwrite(lotusDB.nodup.dt, "Ressources/1_Clean_db_files/LOTU.csv")

# for openbabel
fwrite(lotusDB.nodup.dt[, c("SMILES")], "Ressources/2_Openbabel_conversion/LOTU.smiles", col.names = F)


# COCONUT_DB --------------------------------------------------------------
# Only SMILES
coconutDB <- fread("Ressources/0_Raw_db_files/COCONUT_DB.smi")
coconutDB[, Source := "coconutDB"]
setnames(coconutDB, c("unique_smiles", "coconut_id"), c("SMILES", "ID"))
# no duplicates in coconutDB
fwrite(coconutDB, "Ressources/1_Clean_db_files/COCO.csv")
# for openbabel
fwrite(coconutDB[, c("SMILES")], "Ressources/2_Openbabel_conversion/COCO.smiles", col.names = F)


# CMAUP -------------------------------------------------------------------
# Tab separated and they put tabs in the values of some columns... well done.
# Has InChI, InChIKey, and canonical SMILES
cmaup <- fread("Ressources/0_Raw_db_files/CMAUPv1.0_download_Ingredients_All.txt", sep = "\t")
cmaup[, Source := "CMAUPv1.0"]
setnames(cmaup, c("Ingredient_ID", "V23", "V24", "V25"), c("ID", "InChI", "InChIKey" , "SMILES"))
colnames(cmaup)
cmaup.dt <- cmaup[, c("Source", "ID", "InChI", "InChIKey", "SMILES")]
cmaup.nodup.dt <- cmaup.dt[!duplicated(cmaup.dt[, c("InChI")]), ]

fwrite(cmaup.nodup.dt, "Ressources/1_Clean_db_files/CMAU.csv")
fwrite(cmaup.nodup.dt[, c("InChI")], "Ressources/2_Openbabel_conversion/CMAU.InChI", col.names = F)

# PubChem -----------------------------------------------------------------

# Full PubChem
# Has only SMILES
#pubchem <- fread("/Volumes/4TB/Users/dpflieger/Projects/Gaquerel/Elser/Project_Databases/NP-Databases/Pubchem/CID-InChI-Key.gz")
pubchem_smiles <- fread("Ressources/0_Raw_db_files/CID-SMILES.txt", col.names = c("ID", "SMILES"))

pubchem_tnpa <- fread("Ressources/0_Raw_db_files/PubChem_substance_cache_6q1PEjRIUfRm2tnDW7uQ7YIdQn0MthbHbOINi3fzH4p36iM_summary The Natural Products Atlas.csv")
pubchem_golm <- fread("Ressources/0_Raw_db_files/PubChem_substance_cache_cTbUii3uSFJ_fMBlQh2JS5uo-ciszmogEAVxbAsUY20LDV8_summary-GOLM.csv")
pubchem_bmrb <- fread("Ressources/0_Raw_db_files/PubChem_substance_cache_qO8NU-bIg3S0WgtDiTtCbVCBe-HPnfBjikbrL5FX-S6RTsU_summary-BMRB.csv")
pubchem_metaboworkbench <- fread("Ressources/0_Raw_db_files/PubChem_substance_cache_uP8dRjPKVnZhXFRF1j0daw-Nh-0jGxmaY78C1niuENd4tyw_summaryMetabolomicsWorkbench.csv")

setnames(pubchem_tnpa, c("cid", "sidsrcname"), c("ID", "Source"))
setnames(pubchem_golm, c("cid", "sidsrcname"), c("ID", "Source"))
setnames(pubchem_bmrb, c("cid", "sidsrcname"), c("ID", "Source"))
setnames(pubchem_metaboworkbench, c("cid", "sidsrcname"), c("ID", "Source"))

# Merge with the smiles code from the whole smiles pubchem database
pubchem_tnpa_smiles <- merge(pubchem_tnpa, pubchem_smiles, by = "ID", all.x = TRUE)
pubchem_golm_smiles <- merge(pubchem_golm, pubchem_smiles, by = "ID", all.x = TRUE)
pubchem_bmrb_smiles <- merge(pubchem_bmrb, pubchem_smiles, by = "ID", all.x = TRUE)
pubchem_metaboworkbench_smiles <- merge(pubchem_metaboworkbench, pubchem_smiles, by = "ID", all.x = TRUE)

# LOTS OF MISSING CID... 
pubchem_tnpa_smiles.nodup <- pubchem_tnpa_smiles[!duplicated(pubchem_tnpa_smiles[, c("SMILES")])][!is.na(SMILES)]
pubchem_golm_smiles.nodup <- pubchem_golm_smiles[!duplicated(pubchem_golm_smiles[, c("SMILES")])][!is.na(SMILES)]
pubchem_bmrb_smiles.nodup <- pubchem_bmrb_smiles[!duplicated(pubchem_bmrb_smiles[, c("SMILES")])][!is.na(SMILES)]
pubchem_metaboworkbench_smiles.nodup <- pubchem_metaboworkbench_smiles[!duplicated(pubchem_metaboworkbench_smiles[, c("SMILES")])][!is.na(SMILES)]

fwrite(pubchem_tnpa_smiles.nodup[, c("Source", "ID", "SMILES")], "Ressources/1_Clean_db_files/TNPA.csv")
fwrite(pubchem_golm_smiles.nodup[, c("Source", "ID", "SMILES")], "Ressources/1_Clean_db_files/GOLM.csv")
fwrite(pubchem_bmrb_smiles.nodup[, c("Source", "ID", "SMILES")], "Ressources/1_Clean_db_files/BMRB.csv")
fwrite(pubchem_metaboworkbench_smiles.nodup[, c("Source", "ID", "SMILES")], "Ressources/1_Clean_db_files/MEWO.csv")

# for openbabel
fwrite(pubchem_tnpa_smiles.nodup[, c("SMILES")], "Ressources/2_Openbabel_conversion/TNPA.smiles", col.names = F)
fwrite(pubchem_golm_smiles.nodup[, c("SMILES")], "Ressources/2_Openbabel_conversion/GOLM.smiles", col.names = F)
fwrite(pubchem_bmrb_smiles.nodup[, c("SMILES")], "Ressources/2_Openbabel_conversion/BMRB.smiles", col.names = F)
fwrite(pubchem_metaboworkbench_smiles.nodup[, c("SMILES")], "Ressources/2_Openbabel_conversion/MEWO.smiles", col.names = F)