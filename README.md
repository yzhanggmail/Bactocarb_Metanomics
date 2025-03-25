# Bactocarb_Metanomics R Code #
Code Accompanying paper **"Metabolic plasticity of the gut microbiome in response to diets differing in glycemic load"**.

**Contents:** At current the repository contains the following modules.
   
**Code:**
Code was written in R version 4.3.1. 
1. "Cazyme_MixedEffectModel_testing.R" conducts Mixed Effect Model for the assocaition of CAZymes RNA/DNA ratio in high or low glycemic dietary load (corresponding to Table3 and Figure3B).
2. "MetaPhlanDNA_MixedEffectModel_testing.R" performs Mixed Effect Model for the association of Meta-Phlan DNA data in high or low glycemic dietary load (related to Figure1C). 
3. "Microbiome_Pathway_MixedEffectModel_testing.R" performs Mixed Effect Model for the association of the microbial pathway in either high or low glycemic dietary load (for Table2 and Figure2A), as well as the association with HOMA IR (related to Table4).
4. "Pathway_DNA_RNA_SpearmanCor.R" calculated the weighted Spearman Correlation between metatranscriptomics and metagenomics of pathways prevalent in >80% of participants at the end of HGL and LGL diet periods.

Loading this repository
Please check out this repository using "https://github.com/yzhanggmail/Bactocarb_Metanomics".

git clone --recurse-submodules https://github.com/yzhanggmail/Bactocarb_Metanomics.git
which automatically retrieves all submodules.
