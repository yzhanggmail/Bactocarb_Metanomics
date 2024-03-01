# Bactocarb_Metaomics Data and Code #
Data and Code accompanying "Metabolic plasticity of the gut microbiome in response to diets differing in glycemic load".

Contents
At current the repository contains the following modules.

Data: 
1. "BactoCARB_Clinical_Info_3timepoint.csv" has the clinical data for 80 subjects with 3 timepoints observations including baseline, HGL and LGL samples.
2. "Pathway_organized_filtered_logratio.csv" has the pathway RNA/DNA ratio data in log2 scale. Data was shown present in >80% of observations in the baseline, HGL and LGL metagenomic samples. 147 pathways obtained after removed super-pathway but put back those having diet difference.
3. "Cazymes_organized_filtered_logratio.csv" has CAZymes RNA/DNA ratio data in log2 scale. For RNA/DNA relative expression of bacterial metabolic pathways, we used a relative expression filtering criterion of >80% prevalence, removed all 55 super-pathways, and put back 5 diet-related super-pathways. The CAZyme RNA/DNA ratio data retained 192 variables after we filtered by RNA/DNA ratio>1 of variables identified in 80% of subjects.
4. "MetaPhlanDNA_organized_filtered_log2scale.csv" are the MetaPhlan DNA data in log2 scale. DNA data had 100 species after filtering at a prevalence >80% and 32 species that were <80% prevalence.

Code:
Code was written in R version 4.3.1. 
"Cazyme_MixedEffectModel_testing.R" conduct Mixed Effect Model for the assocaition of CAZymes RNA/DNA ratio in high or low glycemic dietary load (corresponding to Table3 and Figure3B).
"MetaPhlanDNA_MixedEffectModel_testing.R" perform Mixed Effect Model for the association of Meta-Phlan DNA data in high or low glycemic dietary load (related to Figure1C). 
"Microbiome_Pathway_MixedEffectModel_testing.R" perform Mixed Effect Model for the association of the microbial pathway in either high or low glycemic dietary load (for Table2 and Figure2A), as well as the association with HOMA IR (related to Table4).


Loading this repository
Please check out this repository using "https://github.com/yzhanggmail/Bactocarb_Metaomics".

git clone --recurse-submodules https://github.com/yzhanggmail/Bactocarb_Metaomics/edit/main/README.md/Evolution.git
which automatically retrieves all submodules.
