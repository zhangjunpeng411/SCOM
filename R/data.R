#' miRNA expression data  
#' 
#' @docType data 
#' @name ACC_miRNA_Exp
#' @aliases ACC_miRNA_Exp
#' @format ACC_miRNA_Exp: A matrix object with 79 ACC 
#' samples (rows) and 894 miRNAs (columns). 
#' @details The matched ncRNA (including miRNA, lncRNA and pseudogene)
#' and mRNA expression data in Adrenocortical carcinoma (ACC) are obtained 
#' from The Cancer Genome Atlas (TCGA, http://cancergenome.nih.gov/). 
NULL

#' mRNA expression data 
#' 
#' @docType data 
#' @name ACC_mRNA_Exp
#' @aliases ACC_mRNA_Exp
#' @format ACC_mRNA_Exp: A matrix object with 79 ACC 
#' samples (rows) and 19068 mRNAs (columns). 
#' @details The matched ncRNA (including miRNA, lncRNA and pseudogene)
#' and mRNA expression data in Adrenocortical carcinoma (ACC) are obtained 
#' from The Cancer Genome Atlas (TCGA, http://cancergenome.nih.gov/). 
NULL

#' lncRNA expression data 
#' 
#' @docType data 
#' @name ACC_lncRNA_Exp
#' @aliases ACC_lncRNA_Exp
#' @format ACC_lncRNA_Exp: A matrix object with 79 ACC 
#' samples (rows) and 4366 lncRNAs (columns). 
#' @details The matched ncRNA (including miRNA, lncRNA and pseudogene)
#' and mRNA expression data in Adrenocortical carcinoma (ACC) are obtained 
#' from The Cancer Genome Atlas (TCGA, http://cancergenome.nih.gov/). 
NULL

#' pseudogene expression data 
#' 
#' @docType data 
#' @name ACC_pseudogene_Exp
#' @aliases ACC_pseudogene_Exp
#' @format ACC_pseudogene_Exp: A matrix object with 79 ACC 
#' samples (rows) and 11952 pseudogenes (columns). 
#' @details The matched ncRNA (including miRNA, lncRNA and pseudogene)
#' and mRNA expression data in Adrenocortical carcinoma (ACC) are obtained 
#' from The Cancer Genome Atlas (TCGA, http://cancergenome.nih.gov/). 
NULL

#' miRNA-lncRNA and miRNA-mRNA interactions 
#' 
#' @docType data 
#' @name miRTarget_lncR_vs_mR
#' @aliases miRTarget_lncR_vs_mR
#' @format miRTarget_lncR_vs_mR: A data.frame object with
#' 975621 miRNA-lncRNA and miRNA-mRNA interactions.
#' @details The putative miRNA-mRNA interactions (762,540) are acquired 
#' from miRTarBase v9.0 and TarBase v8.0, and the putative miRNA-lncRNA 
#' interactions (225,063) are obtained from LncBase v2.0 and NPInter 4.0.
#' @references Huang H, Lin Y-C, Cui S, Huang Y, Tang Y, Xu J, et al. 
#' miRTarBase update 2022: an informative resource for experimentally 
#' validated miRNA-target interactions. Nucleic Acids Res. 
#' 2022;50: D222–D230. doi:10.1093/nar/gkab1079
#' @references Karagkouni D, Paraskevopoulou MD, Chatzopoulos S, Vlachos IS, 
#' Tastsoglou S, Kanellos I, et al. DIANA-TarBase v8: a decade-long 
#' collection of experimentally supported miRNA-gene interactions. 
#' Nucleic Acids Res. 2018;46: D239–D245. doi:10.1093/nar/gkx1141
#' @references Paraskevopoulou MD, Vlachos IS, Karagkouni D, Georgakilas G, 
#' Kanellos I, Vergoulis T, et al. DIANA-LncBase v2: indexing microRNA 
#' targets on non-coding transcripts. Nucleic Acids Res. 2016;44: D231-238. 
#' doi:10.1093/nar/gkv1270
#' @references Teng X, Chen X, Xue H, Tang Y, Zhang P, Kang Q, et al. 
#' NPInter v4.0: an integrated database of ncRNA interactions. Nucleic Acids Res. 
#' 2020;48: D160–D165. doi:10.1093/nar/gkz969
NULL

#' miRNA-pseudogene and miRNA-mRNA interactions 
#' 
#' @docType data 
#' @name miRTarget_pseudo_vs_mR
#' @aliases miRTarget_pseudo_vs_mR
#' @format miRTarget_pseudo_vs_mR: A data.frame object with
#' 854141 miRNA-pseudogene and miRNA-mRNA interactions.
#' @details The putative miRNA-mRNA interactions (762,540) are acquired 
#' from miRTarBase v9.0 and TarBase v8.0, and the putative miRNA-pseudogene 
#' interactions (91,619) are from ENCORI (the previous version is starBase).
#' @references Huang H, Lin Y-C, Cui S, Huang Y, Tang Y, Xu J, et al. 
#' miRTarBase update 2022: an informative resource for experimentally 
#' validated miRNA-target interactions. Nucleic Acids Res. 
#' 2022;50: D222–D230. doi:10.1093/nar/gkab1079
#' @references Karagkouni D, Paraskevopoulou MD, Chatzopoulos S, Vlachos IS, 
#' Tastsoglou S, Kanellos I, et al. DIANA-TarBase v8: a decade-long 
#' collection of experimentally supported miRNA-gene interactions. 
#' Nucleic Acids Res. 2018;46: D239–D245. doi:10.1093/nar/gkx1141
#' @references Li JH, Liu S, Zhou H, Qu LH, Yang JH. 
#' starBase v2.0: decoding miRNA-ceRNA, miRNA-ncRNA and protein-RNA 
#' interaction networks from large-scale CLIP-Seq data. Nucleic Acids Res. 
#' 2014;42: D92-97. doi:10.1093/nar/gkt1248
NULL