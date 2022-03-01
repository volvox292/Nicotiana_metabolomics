library(tidyverse)
library(data.table)
library(plyr)

converted_dir = "Ressources/2_Openbabel_conversion/converted"

# SUBS --------------------------------------------------------------------
SUBS <- fread("Ressources/1_Clean_db_files/SUBS.csv")
SUBS.inchi <- fread(file.path(converted_dir, "SUBS.inchi"), sep = "\t", header = F, col.names = c("InChI"))
SUBS.inchikey <- fread(file.path(converted_dir, "SUBS.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
SUBS.dt <- bind_cols(SUBS, SUBS.inchi, SUBS.inchikey)

# MENP --------------------------------------------------------------------
MENP <- fread("Ressources/1_Clean_db_files/MENP.csv")
MENP.inchi <- fread(file.path(converted_dir, "MENP.inchi"), sep = "\t", header = F, col.names = c("InChI"))
MENP.inchikey <- fread(file.path(converted_dir, "MENP.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
MENP.dt <- bind_cols(MENP, MENP.inchi, MENP.inchikey)

# LOTU --------------------------------------------------------------------
LOTU <- fread("Ressources/1_Clean_db_files/LOTU.csv")
LOTU.inchi <- fread(file.path(converted_dir, "LOTU.inchi"), sep = "\t", header = F, col.names = c("InChI"))
LOTU.inchikey <- fread(file.path(converted_dir, "LOTU.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
LOTU.dt <- bind_cols(LOTU, LOTU.inchi, LOTU.inchikey)

# COCO --------------------------------------------------------------------
COCO <- fread("Ressources/1_Clean_db_files/COCO.csv")
COCO.inchi <- fread(file.path(converted_dir, "COCO.inchi"), sep = "\t", header = F, col.names = c("InChI"))
COCO.inchikey <- fread(file.path(converted_dir, "COCO.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
COCO.dt <- bind_cols(COCO, COCO.inchi, COCO.inchikey)

# GOLM --------------------------------------------------------------------
GOLM <- fread("Ressources/1_Clean_db_files/GOLM.csv")
GOLM.inchi <- fread(file.path(converted_dir, "GOLM.inchi"), sep = "\t", header = F, col.names = c("InChI"))
GOLM.inchikey <- fread(file.path(converted_dir, "GOLM.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
GOLM.dt <- bind_cols(GOLM, GOLM.inchi, GOLM.inchikey)

# MEWO --------------------------------------------------------------------
MEWO <- fread("Ressources/1_Clean_db_files/MEWO.csv")
MEWO.inchi <- fread(file.path(converted_dir, "MEWO.inchi"), sep = "\t", header = F, col.names = c("InChI"))
MEWO.inchikey <- fread(file.path(converted_dir, "MEWO.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
MEWO.dt <- bind_cols(MEWO, MEWO.inchi, MEWO.inchikey)

# TNPA --------------------------------------------------------------------
TNPA <- fread("Ressources/1_Clean_db_files/TNPA.csv")
TNPA.inchi <- fread(file.path(converted_dir, "TNPA.inchi"), sep = "\t", header = F, col.names = c("InChI"))
TNPA.inchikey <- fread(file.path(converted_dir, "TNPA.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
TNPA.dt <- bind_cols(TNPA, TNPA.inchi, TNPA.inchikey)

# BMRB --------------------------------------------------------------------
BMRB <- fread("Ressources/1_Clean_db_files/BMRB.csv")
BMRB.inchi <- fread(file.path(converted_dir, "BMRB.inchi"), sep = "\t", header = F, col.names = c("InChI"))
BMRB.inchikey <- fread(file.path(converted_dir, "BMRB.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
BMRB.dt <- bind_cols(BMRB, BMRB.inchi, BMRB.inchikey)

# KNAP --------------------------------------------------------------------
KNAP <- fread("Ressources/1_Clean_db_files/KNAP.csv")
# openbabel when converting to smiles makes weird file, we have two columns
KNAP.smiles <- fread(file.path(converted_dir, "KNAP.smiles"), sep = "\t", header = F, col.names = c("SMILES", "nothing"))
KNAP.smiles[, nothing := NULL]
KNAP.inchikey <- fread(file.path(converted_dir, "KNAP.inchikey"), sep = "\t", header = F, col.names = c("InChIKey"))
KNAP.dt <- bind_cols(KNAP, KNAP.smiles, KNAP.inchikey)

UNPD.dt <- fread("Ressources/1_Clean_db_files/UNPD.csv")
CMAU.dt <- fread("Ressources/1_Clean_db_files/CMAU.csv")

# Merge databanks ---------------------------------------------------------
dt.list <- list(SUBS.dt, UNPD.dt, KNAP.dt, MENP.dt, COCO.dt, LOTU.dt, CMAU.dt, TNPA.dt, GOLM.dt, BMRB.dt, MEWO.dt)

lapply(dt.list, colnames)

total.dt <- rbindlist(dt.list, use.names = T)

# Number of compounds in our database 
nrow(total.dt)

# Now we need to get missing SMILES, InChI and InChIKey
# Then reduce again for redundant InChI
total.nodup.inchi.dt <- total.dt[!duplicated(total.dt, by = "InChI")]
nrow(total.nodup.inchi.dt)

# Same InChI but not same SMILES
total.dt[InChI == "InChI=1S/C15H14O9/c1-8(16)23-12(14(19)20)13(15(21)22)24-11(18)7-4-9-2-5-10(17)6-3-9/h2-7,12-13,17H,1H3,(H,19,20)(H,21,22)/b7-4+/t12-,13-/m1/s1"]

# Duplicated ?
# SMILES 133228
# INCHI 283547

# SQLite DB ---------------------------------------------------------------
library(DBI)

# Create a "Database" folder for the sqlite database
dir.create("Database/", showWarnings = FALSE)

# Database name 
SQLITE_DB = "compounds.sqlite"

# Init the database
mydb <- dbConnect(RSQLite::SQLite(), file.path("Database", SQLITE_DB))
dbListTables(mydb)

# Load data into the database in a "compounds" table
dbWriteTable(mydb, "compounds", total.nodup.inchi.dt)

# Query and extract as csv from the database
my.data.frame <- as.data.table(dbGetQuery(mydb, "SELECT ID, InchI FROM compounds"))
fwrite(my.data.frame, "nr_compounds_db.txt", sep = " ", col.names = F)


#### Added the 25.08.2021 ####
# TESTING - NEW COMPOUNDS FROM A NICOTIANA REVIEW
nicotiana_review <- fread("/Volumes/4TB/Users/dpflieger/Projects/Gaquerel/Elser/Compounds_databank/Ressources/tmp/jassbi_structures.txt", header = F)
setnames(nicotiana_review, c("ID", "InChI"))
dim(nicotiana_review)


# Compounds not matching exactly to what we have in the database
not_matching_exact <- nicotiana_review[!InChI %in% my.data.frame$InChI]
not_matching_exact[, short_InChI := substr(InChI, 0, 50)]

my.data.frame[, short_InChI := substr(InChI, 0, 50)]

res <- merge(not_matching_exact, my.data.frame, by = "short_InChI", all.x = T)

fwrite(res, "/Volumes/4TB/Users/dpflieger/tmp/tmp.csv", sep = "\t")

unique(res[is.na(ID.y)]$short_InChI)

       