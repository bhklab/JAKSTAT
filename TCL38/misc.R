dir <- "~/Documents/JAKSTAT/tcl38/"

library(dplyr) 
library(tidyr)
library(readxl)
library(stringr)
library(ghql)
library(jsonlite)
library(rlist)
library(PharmacoGx)
library(plotly)

source(file='tcl38_functions.R')

pset.list <- readRDS("./misc/data/tcl38_pset_list.rds")

cell_drug_combinations <- get_cell_drug_combinations(pset.list)

cells <- levels(cell_drug_combinations$cellid)
drugs <- levels(cell_drug_combinations$drugid)

cells_df <- data.frame(matrix(data=NA, ncol=0, nrow=length(cells)))
cells_df["cellid"] <- cells
for(dataset.name in names(pset.list)){
  dataset_cells <- pset.list[[dataset.name]]@cell
  cells_df[dataset.name] <- ifelse(
    is.na(match(
      cells_df$cellid,
      dataset_cells$cellid,
    )),
    NA, "Yes"
  )
}

drugs_df <- data.frame(matrix(data=NA, ncol=0, nrow=length(drugs)))
drugs_df["drugid"] <- drugs
for(dataset.name in names(pset.list)){
  dataset_cells <- pset.list[[dataset.name]]@drug
  drugs_df[dataset.name] <- ifelse(
    is.na(match(
      drugs_df$drugid,
      dataset_cells$drugid,
    )),
    NA, "Yes"
  )
}

mol_profs <- c()
for(dataset.name in names(pset.list)){
  mol_profs <- c(mol_profs, names(pset.list[[dataset.name]]@molecularProfiles))
}
mol_profs <- unique(mol_profs)

mol_profs_df <- data.frame(matrix(data=NA, ncol=0, nrow=length(mol_profs)))
mol_profs_df["molecular.profile"] <- mol_profs
for(dataset.name in names(pset.list)){
  dataset_mol_prof <- names(pset.list[[dataset.name]]@molecularProfiles)
  mol_profs_df[dataset.name] <- ifelse(
    is.na(match(
      mol_profs,
      dataset_mol_prof,
    )),
    NA, "Yes"
  )
}
mol_profs_df$molecular.profile[mol_profs_df$molecular.profile == "rna"] <- "rna(microarray)"

summary <- cells_df
summary <- summary[order(summary$cellid), ]
tmp <- summary[, names(pset.list)]
summary["datasets"] <- apply(tmp, 1, function(row){paste(colnames(tmp)[ !is.na(row) ], collapse = ", ")})
summary[names(pset.list)] <- NULL

tmp <- cell_drug_combinations[rowSums(is.na(cell_drug_combinations)) <= 4, ]
tmp <- data.frame(table(tmp$cellid[tmp$cellid %in% summary$cellid]))
tmp <- tmp[order(tmp$Var1), ]
summary["num.drugs.tested"] <- tmp$Freq
summary[mol_profs[order(mol_profs)]] <- NA

for(dataset.name in names(pset.list)){
  for(molprof.name in names(pset.list[[dataset.name]]@molecularProfiles)){
    tmp <- data.frame(pset.list[[dataset.name]]@molecularProfiles[[molprof.name]]@colData)
    origin <- summary[[molprof.name]]
    incoming <- ifelse(summary$cellid %in% tmp$cellid, "Yes", NA)
    summary[[molprof.name]] <- dplyr::coalesce(origin, incoming)
  }
}

colnames(summary)[colnames(summary) == "rna"] <- "rna(microarray)"

write.csv(summary, file=paste0(dir, "cell_molecular_data_summary.csv"))

# Assign MoA to each drug
chembl <- read.csv(paste0(dir, "data-moa-filtered.csv"))
chembl <- chembl[, c("compound_name", "mechanism_of_action", "action_type", "standard_inchi", "standard_inchi_key", "canonical_smiles")]
pharmacodb.drugs <- read.csv(paste0(dir, "drugs_pharmacodb.csv"))
tcl38.drugs <- read.csv(paste0(dir, "drugs_summary.csv"))
tcl38.drugs[c("X", "CCLE", "GDSC1", "GDSC2", "gCSI", "CTRPv2")] <- NULL
tcl38.drugs$drugid <- tcl38.drugs$drugid[order(tcl38.drugs$drugid)]
pharmacodb.drugs <- pharmacodb.drugs[tolower(pharmacodb.drugs$name) %in% tolower(tcl38.drugs$drugid), ]
colnames(pharmacodb.drugs)[colnames(pharmacodb.drugs) == "name"] <- "drugid"
tcl38.drugs$drugid[!(tcl38.drugs$drugid %in% pharmacodb.drugs$name)] #11 of them that could not be found

tcl38.drugs <- merge(x=tcl38.drugs, y=pharmacodb.drugs[, c("drugid", "smiles", "inchikey")], by="drugid", all.x=TRUE)

inchikeys <- tcl38.drugs$inchikey[(!is.na(tcl38.drugs$inchikey) & tcl38.drugs$inchikey != "NULL" & tcl38.drugs$inchikey != "None")]
chembl.filtered <- chembl[chembl$standard_inchi_key %in% inchikeys, ]
colnames(chembl.filtered)[colnames(chembl.filtered) == "standard_inchi_key"] <- "inchikey"
tcl38.drugs <- merge(x=tcl38.drugs, y=chembl.filtered[, c("mechanism_of_action", "action_type", "inchikey", "canonical_smiles")], by="inchikey", all.x=TRUE)
tcl38.drugs <- tcl38.drugs[, c("drugid", "mechanism_of_action", "action_type", "inchikey", "smiles", "canonical_smiles")]
tcl38.drugs$action_type[is.na(tcl38.drugs$action_type)] <- "Unidentified"
tcl38.drugs$mechanism_of_action[is.na(tcl38.drugs$mechanism_of_action)] <- "N/A"

tcl38.drugs.subsetted <- tcl38.drugs[, c("drugid", "mechanism_of_action", "action_type")]
tcl38.drugs.subsetted <- tcl38.drugs.subsetted[!duplicated(tcl38.drugs.subsetted[c( "drugid", "mechanism_of_action", "action_type")]), ]
tcl38.drugs.subsetted <- tcl38.drugs.subsetted[tcl38.drugs.subsetted$mechanism_of_action != "N/A", ]


# sunburst diagram
moas.summary <- read.csv(paste0(dir, "moas_summary.csv"))
sunburst.data <- merge(x=tcl38.drugs.subsetted, y=moas.summary[, c("mechanism_of_action", "inhibitor_subgroup")], by="mechanism_of_action", all=TRUE)

action.types <- data.frame(table(sunburst.data$action_type))
inhibitor.sub <- data.frame(table(sunburst.data$inhibitor_subgroup))
inhibitor.sub <- inhibitor.sub[inhibitor.sub$Var1 != "", ]
moas <- data.frame(table(sunburst.data$mechanism_of_action))

labels <- c(levels(action.types$Var1), levels(inhibitor.sub$Var1), levels(moas$Var1))
labels <- labels[labels != ""]
values <- c(action.types$Freq, inhibitor.sub$Freq, moas$Freq)

labels <- c(paste(labels, "<br>(", values, ")"), "Mechanisms<br> of<br> Action")
values <- c(values, sum(action.types$Freq))

parents <- c(as.vector(matrix("Mechanisms<br> of<br> Action", nrow=length(action.types$Var1))))
parents <- c(parents, c(as.vector(matrix("INHIBITOR <br>( 163 )", nrow=length(inhibitor.sub$Var1)))))

moas.action.type <- sunburst.data[, c("mechanism_of_action", "inhibitor_subgroup", "action_type")]
moas.action.type <- moas.action.type[!duplicated(moas.action.type[c("mechanism_of_action", "inhibitor_subgroup", "action_type")]), ]
action.type.parents <- c()
for(moa in moas$Var1){
  inhibitor.subgroup <- moas.action.type[moas.action.type$mechanism_of_action == moa, ]$inhibitor_subgroup
  if(inhibitor.subgroup != ""){
    num <- inhibitor.sub[inhibitor.sub$Var1 == inhibitor.subgroup, ]$Freq
    action.type.parents <- c(action.type.parents, paste(inhibitor.subgroup, "<br>(", num, ")"))
  }else{
    action.type <- moas.action.type[moas.action.type$mechanism_of_action == moa, ]$action_type
    num <- action.types[action.types$Var1 == action.type, ]$Freq
    # action.type.parents <- c(action.type.parents, moas.action.type[moas.action.type$mechanism_of_action == moa, ]$action_type)
    action.type.parents <- c(action.type.parents, paste(action.type, "<br>(", num, ")"))
  }
}
parents <- c(parents, action.type.parents, "")

fig <- plot_ly(
  labels = labels,
  parents = parents,
  values = values,
  type = 'sunburst',
  branchvalues = 'total',
  insidetextorientation='radial'
)
