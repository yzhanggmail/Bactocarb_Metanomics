# Bactocarb_Metaomics Data and Code #
Data and Code Accompanying "Metabolic plasticity of the gut microbiome in response to diets differing in glycemic load".

**Contents:** At current the repository contains the following modules.

Data: 
1. "BactoCARB_Clinical_Info_3timepoint.csv" contains clinical data for 80 subjects with observations at three timepoints including baseline, HGL and LGL samples.
2. "Pathway_organized_filtered_logratio.csv" includes pathway RNA/DNA ratio data for HGL and LGL observations in log2 scale. The data were present in >80% of observations in the baseline, HGL and LGL metagenomic samples. A total of 147 pathways werer obtained after removing super-pathway but those showing diet difference were included.
3. "Cazymes_organized_filtered_logratio.csv" contains CAZymes RNA/DNA ratio data for HGL and LGL observations in log2 scale. For RNA/DNA relative expression of bacterial metabolic pathways, a relative expression filtering criterion of >80% prevalence was used. All 55 super-pathways were removed, and 5 diet-related super-pathways were reinstated. The CAZyme RNA/DNA ratio data retained 192 variables after filtering by RNA/DNA ratio>1 of variables identified in 80% of subjects.
4. "MetaPhlanDNA_organized_filtered_log2scale.csv" includes MetaPhlan DNA data for HGL and LGL observations in log2 scale. The DNA data consisted of 100 species after filtering at a prevalence >80%.
Code:
Code was written in R version 4.3.1. 
1. "Cazyme_MixedEffectModel_testing.R" conduct Mixed Effect Model for the assocaition of CAZymes RNA/DNA ratio in high or low glycemic dietary load (corresponding to Table3 and Figure3B).
2. "MetaPhlanDNA_MixedEffectModel_testing.R" perform Mixed Effect Model for the association of Meta-Phlan DNA data in high or low glycemic dietary load (related to Figure1C). 
3. "Microbiome_Pathway_MixedEffectModel_testing.R" perform Mixed Effect Model for the association of the microbial pathway in either high or low glycemic dietary load (for Table2 and Figure2A), as well as the association with HOMA IR (related to Table4).


Loading this repository
Please check out this repository using "https://github.com/yzhanggmail/Bactocarb_Metaomics".

git clone --recurse-submodules https://github.com/yzhanggmail/Bactocarb_Metaomics/edit/main/README.md/Evolution.git
which automatically retrieves all submodules.
