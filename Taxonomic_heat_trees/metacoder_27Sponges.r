Rstudio (R 4.4.2, RStudio version 2024.12.1+563)

library("ggplot2")
library("dplyr")
library("phyloseq")
library("metacoder")
library("readr")

library("readxl")
library("tibble")


# load sheets from excel file. 
# Note that the taxonomy table has NA where taxon was unassigned (eg s__ or g__) and placeholder names starting with a number have been moved to Xplaceholdername. 
# Note that the OTU_table_relabund has been converted to 100 000 000 reads per sponge from the % relative abundances calculated by BioSAK magabund.

otu_nano <- read_excel("metacoder_inputs.xlsx", sheet="OTU_table_relabund")

tax_nano <- read_excel("metacoder_inputs.xlsx", sheet="taxonomy_table")

samples_nano <- read_excel("metacoder_inputs.xlsx", sheet="metadata_table")


# assign the row names
otu_nano <- otu_nano %>% tibble::column_to_rownames("ASV")
tax_nano <- tax_nano %>% tibble::column_to_rownames("ASV")
samples_nano <- samples_nano %>% tibble::column_to_rownames("Samples")

# convert to matrix (not needed for metadata samples file)
otu_nano <- as.matrix(otu_nano)
tax_nano <- as.matrix(tax_nano)

# Making input for phyloseq
OTU = otu_table(otu_nano, taxa_are_rows = TRUE)
TAX = tax_table(tax_nano)
samples = sample_data(samples_nano)

# Make phyloseq object
sponge27taxa <- phyloseq(OTU, TAX, samples) 		
sponge27taxa	

#	phyloseq-class experiment-level object
#	otu_table()   OTU Table:         [ 302 taxa and 26 samples ]
#	sample_data() Sample Data:       [ 26 samples by 35 sample variables ]
#	tax_table()   Taxonomy Table:    [ 302 taxa by 8 taxonomic ranks ] 		


obj <- parse_phyloseq(sponge27taxa)

metadata_obj <- obj$data$sample_data

obj$data$rel_abd <- calc_obs_props(obj, "otu_table", other_cols = T) # Convert to rel abundance

obj$data$tax_rel_abd <- calc_taxon_abund(obj, "rel_abd")  # Now calculate rel abundance per taxon



obj$data$obs_prop <- calc_obs_props(obj, "otu_table", other_cols = T, out_names = c("Aphro05obs", "Aphro18obs", "Boloobs", "Calyxobs", "Charaobs", "Chondobs", "Farreaobs", "Ggar08obs", "Ggar09obs", "Gmega16obs", "Gmega17obs", "GparvCASobs", "GparvKARobs", "GparvNORobs", "Homaxobs", "Paratiobs", "Pcalic07obs", "Pcalic13obs", "Pcalic14obs", "Phironobs", "Phycoobs", "Povistobs", "PventBARobs", "PventNORobs", "StryBARobs", "StryCANobs", "StryNORobs"))


obj$data$taxrelabd <- calc_obs_props(obj, "tax_rel_abd", other_cols = T, out_names = c("Aphro05taxrelabd", "Aphro18taxrelabd", "Bolotaxrelabd", "Calyxtaxrelabd", "Charataxrelabd", "Chondtaxrelabd", "Farreataxrelabd", "Ggar08taxrelabd", "Ggar09taxrelabd", "Gmega16taxrelabd", "Gmega17taxrelabd", "GparvCAStaxrelabd", "GparvKARtaxrelabd", "GparvNORtaxrelabd", "Homaxtaxrelabd", "Paratitaxrelabd", "Pcalic07taxrelabd", "Pcalic13taxrelabd", "Pcalic14taxrelabd", "Phirontaxrelabd", "Phycotaxrelabd", "Povisttaxrelabd", "PventBARtaxrelabd", "PventNORtaxrelabd", "StryBARtaxrelabd", "StryCANtaxrelabd", "StryNORtaxrelabd"))

obj$data$tax_rel_abd <- calc_taxon_abund(obj, "rel_abd", out_names = c("Aphro05tax_rel_abd", "Aphro18tax_rel_abd", "Bolotax_rel_abd", "Calyxtax_rel_abd", "Charatax_rel_abd", "Chondtax_rel_abd", "Farreatax_rel_abd", "Ggar08tax_rel_abd", "Ggar09tax_rel_abd", "Gmega16tax_rel_abd", "Gmega17tax_rel_abd", "GparvCAStax_rel_abd", "GparvKARtax_rel_abd", "GparvNORtax_rel_abd", "Homaxtax_rel_abd", "Paratitax_rel_abd", "Pcalic07tax_rel_abd", "Pcalic13tax_rel_abd", "Pcalic14tax_rel_abd", "Phirontax_rel_abd", "Phycotax_rel_abd", "Povisttax_rel_abd", "PventBARtax_rel_abd", "PventNORtax_rel_abd", "StryBARtax_rel_abd", "StryCANtax_rel_abd", "StryNORtax_rel_abd"))


taxonFulllist <- obj$taxon_names()
transposed <- t(taxonFulllist)
transposed2 <- t(transposed)


### For each SPONGE (where SPONGE = Aphro05, Aphro18, Bolo, Calyx, Chara, Chond, Farrea, Ggar08, Ggar09, Gmega16, Gmega17, GparvCAS, GparvKAR, GparvNOR, Homax, Parati, Pcalic07, Pcalic13, Pcalic14, Phiron, Phyco, Povist, PventBAR, PventNOR, StryBAR, StryCAN, StryNOR) use the following commands to produce the plot

nameslist_SPONGE <- filter(obj$data$tax_rel_abd, SPONGEtax_rel_abd == "0")
nameslist_SPONGE_ids <-nameslist_SPONGE$taxon_id
matching_rows_SPONGE <- transposed2[row.names(transposed2) %in% nameslist_SPONGE_ids, ] # Prep so that only taxa present in that sponge will be labelled

set.seed(5)
obj %>% filter_taxa(taxon_ranks == "Genus", supertaxa = TRUE) %>% heat_tree(node_label = ifelse(taxon_names %in% matching_rows_SPONGE, "", taxon_names), node_size = SPONGEtax_rel_abd, node_color = SPONGEtax_rel_abd, edge_color = SPONGEtax_rel_abd, initial_layout = "fr", layout = "da", node_color_interval = c(0, 1), node_size_range = c(0.005, 0.05), node_label_size_range = c(0.008, 0.04), output_file = "plot_SPONGE_scaled_seed5.pdf")




# for example

nameslist_Chond <- filter(obj$data$tax_rel_abd, Chondtax_rel_abd == "0")
nameslist_Chond_ids <-nameslist_Chond$taxon_id
matching_rows_Chond <- transposed2[row.names(transposed2) %in% nameslist_Chond_ids, ]

set.seed(5)
obj %>% filter_taxa(taxon_ranks == "Genus", supertaxa = TRUE) %>% heat_tree(node_label = ifelse(taxon_names %in% matching_rows_Chond, "", taxon_names), node_size = Chondtax_rel_abd, node_color = Chondtax_rel_abd, edge_color = Chondtax_rel_abd, initial_layout = "fr", layout = "da", node_color_interval = c(0, 1), node_size_range = c(0.005, 0.05), node_label_size_range = c(0.008, 0.04), output_file = "plot_Chond_scaled_seed5.pdf")
