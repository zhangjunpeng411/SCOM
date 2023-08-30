# :hammer: SCOM
**Pan-cancer characterization of ncRNA synergistic competition uncovers potential carcinogenic biomarkers**

# Installation
```{r echo=FALSE, results='hide', message=FALSE}
install.packages("devtools")
library(devtools)
install_github("zhangjunpeng411/SCOM")
```

## :boom: Background
The ceRNA regulation of malignant tumors involves the interactions between multiple non-coding and coding genes. In the ceRNA network in which ncRNA is involved, the competitive relationships between ncRNAs and target genes are not one-to-one, but many-to-many. Furthermore, the entire cell system presents the characteristics of functional modularity, and each gene functional module (that is, subsets of genes are more frequently interact or crosstalk with each other than estimated by chance) is responsible for specific biological functions in human complex diseases, including malignant tumors. These findings indicate that ncRNAs acting as ceRNAs not independently, but synergistically play important roles in malignant tumors. 

In this work, to uncover potential carcinogenic biomarkers, we present a **S**ynergistic **COM**petition (**SCOM**) framework to characterize ncRNA synergistic competition in Pan-cancer. Based on the proposed ncRNA synergistic competition hypothesis, **SCOM** applies three criteria (**significantly synergistic competition for mRNAs**, **significantly positive correlation**, and **significantly sensitive correlation conditioning on synergistically competed mRNAs**) to predict ncRNA synergistic competition network from gene (ncRNAs acting as ceRNAs, and target mRNAs) expression data and predicted ncRNA-related ceRNA networks. Driven by Pan-cancer transcriptomics data and putative miRNA-target interactions, **SCOM** lays a foundation for establishing the relationship between ncRNA and its synergistic competition and malignant tumors, and can help to elucidate ncRNA synergistic competition mechanisms in malignant tumors.

A schematic illustration of **SCOM** is shown in the folowing.

<p align="center">
  <img src="https://github.com/zhangjunpeng411/SCOM/blob/main/SCOM_schematic_illustration.png" alt="SCOM schematic illustration" border="0.1">
</p>

**a.** An illustration of the ncRNA synergistic competition hypothesis. 

**b.** For each tumor type in Pan-cancer, **SCOM** firstly predicts ncRNA-related ceRNA network by incorporating gene expression data in matched tumor tissues and priori information of miRNA-target interactions. By integrating gene expression data in matched tumor tissues and predicted ncRNA-related ceRNA networks, **SCOM** further infers ncRNA synergistic competition network for each tumor type in Pan-cancer.

## :zap: Quick example to use SCOM
For inferring ncRNA synergistic competition, users should prepare matched miRNA, lncRNA, pseudogene and mRNA expression data and putative miRNA-target interactions. Users can use the following scripts to infer ncRNA synergistic competition. 

```{r echo=FALSE, results='hide', message=FALSE}
## Load SCOM package
library(SCOM)

## Load prepared datasets in SCOM package
data(ACC)

## Identifying ceRNA interactions and ncRNA synergistic competitions in ACC samples
library(doParallel)
library(SPONGE)
num.cores <- 2
cores <- makeCluster(num.cores)
registerDoParallel(cores)

pre_null_model_ACC <- sponge_build_null_model(number_of_datasets = 1e+03, 
	                                  number_of_samples = nrow(ACC_miRNA_Exp))
ceRNet_lncR_vs_mR_ACC <- SC(ACC_miRNA_Exp, ACC_lncRNA_Exp, ACC_mRNA_Exp, miRTarget_lncR_vs_mR, null_model = pre_null_model_ACC)
ceRNet_pseudo_vs_mR_ACC <- SC(ACC_miRNA_Exp, ACC_pseudogene_Exp, ACC_mRNA_Exp, miRTarget_pseudo_vs_mR, null_model = pre_null_model_ACC)
ceRNet_ACC <- rbind(ceRNet_lncR_vs_mR_ACC, ceRNet_pseudo_vs_mR_ACC)
SCOMNet_ACC <- SCOM(ceRNet_ACC, cbind(ACC_lncRNA_Exp, ACC_pseudogene_Exp), ACC_mRNA_Exp, null_model = pre_null_model_ACC)
```    

## License
GPL-3
