source("omega/OmegaMain.R")
source("pathways/HeatmapFunction.R")
source("pathways/KEGGtable.R")


kegg.pathway <- function(pathway.description) {
  orthogroups <- ko.table %>%
    filter(pathway == pathway.description | category2 == pathway.description | category1 == pathway.description) %>%
    pull(new_orthogroup) %>%
    strsplit(", ") %>%
    unlist() %>%
    unique()
  
  png(paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/", sample.size, "/output/KEGG_pathways/", pathway.description, ".png") , width = 1200, height = 1600)
  gc.heatmap(orthogroups, pathway.description)
  dev.off() 
  
  pathway <- orthogroup.sample %>%
    filter(new_orthogroup %in% orthogroups) %>%
    mutate(new_orthogroup = factor(new_orthogroup, levels = orthogroups)) %>%
    arrange(new_orthogroup) 
  
  pathway$new_orthogroup <- as.numeric(as.character(pathway$new_orthogroup))
  
  pathway.omega <<- pathway %>%
    left_join(M0.omega[c("new_orthogroup", "omega", "number_family_with_loss")], by = "new_orthogroup") %>%
    left_join(dataset.kegg[c("new_orthogroup", "EC")], by = "new_orthogroup") %>%
    select(omega, EC, n, number_family_with_loss, everything()) 
}



#carbohydrate metabolism
kegg.pathway("Glycolysis_or_Gluconeogenesis")
kegg.pathway("Citrate_cycle_TCA_cycle")
kegg.pathway("Pentose_phosphate_pathway") #M00007 non-oxidative phase
kegg.pathway("Pyruvate_metabolism")
kegg.pathway("Glyoxylate_and_dicarboxylate_metabolism") #M00012 Glyoxylate cycle
kegg.pathway("C5_Branched_dibasic_acid_metabolism")

#energy metabolism
kegg.pathway("Oxidative_phosphorylation")
kegg.pathway("Nitrogen_metabolism")
kegg.pathway("Sulfur_metabolism")

#lipid metabolism
kegg.pathway("Fatty_acid_biosynthesis")
kegg.pathway("Glycerolipid_metabolism")
kegg.pathway("Glycerophospholipid_metabolism")

#nucleotide metabolism
kegg.pathway("Purine_metabolism")
kegg.pathway("Pyrimidine_metabolism")

#amino acid metabolism
#kegg.pathway("Biosynthesis_of_amino_acids")
kegg.pathway("Alanine_aspartate_and_glutamate_metabolism")
kegg.pathway("Glycine_serine_and_threonine_metabolism")
kegg.pathway("Cysteine_and_methionine_metabolism")
kegg.pathway("Valine_leucine_and_isoleucine_biosynthesis")
kegg.pathway("Lysine_biosynthesis")
kegg.pathway("Arginine_biosynthesis")
kegg.pathway("Histidine_metabolism")
kegg.pathway("Phenylalanine_tyrosine_and_tryptophan_biosynthesis")

#glycan biosynthesis and metabolism
kegg.pathway("Peptidoglycan_biosynthesis")


#Metabolism of cofactors and vitamins
#kegg.pathway("Metabolism of cofactors and vitamins") #category2
kegg.pathway("Thiamine_metabolism") #B1
kegg.pathway("Riboflavin_metabolism") #B2
kegg.pathway("Nicotinate_and_nicotinamide_metabolism") #B3
kegg.pathway("Pantothenate_and_CoA_biosynthesis") #B5
kegg.pathway("Vitamin_B6_metabolism") #B6
kegg.pathway("Biotin_metabolism") #B7
kegg.pathway("Folate_biosynthesis") #B9
kegg.pathway("One_carbon_pool_by_folate")
kegg.pathway("Lipoic_acid_metabolism")
kegg.pathway("Porphyrin_metabolism") #M00846 siroheme biosynthesis
kegg.pathway("Ubiquinone_and_other_terpenoid_quinone_biosynthesis")

#Metabolism of terpenoids and polyketides
kegg.pathway("Terpenoid_backbone_biosynthesis")

#Transcription and Translation
kegg.pathway("RNA_polymerase")
kegg.pathway("Ribosome")
kegg.pathway("Aminoacyl_tRNA_biosynthesis") #loss of glnS in some strains

#Folding sorting and degradation
kegg.pathway("Protein_export")
kegg.pathway("RNA_degradation")

#Replication and repair
kegg.pathway("Replication and repair") #category2
kegg.pathway("DNA_replication")
kegg.pathway("Base_excision_repair")
kegg.pathway("Mismatch_repair")
kegg.pathway("Homologous_recombination")

#Membrane transport
kegg.pathway("ABC_transporters")
kegg.pathway("Bacterial_secretion_system")

#Signal transduction
kegg.pathway("Two_component_system")

#Cell growth and death
kegg.pathway("Cell growth and death") #category2

