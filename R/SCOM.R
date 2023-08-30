#' Identifying ceRNA interactions using Sensitivity Correlation (SC) method and null model.
#' 
#' @title SC
#' @param miRExp miRNA expression data, rows are samples, columns are miRNAs.
#' @param ceRExp ceRNA (lncRNAs or pseudogenes) expression data, rows are samples, columns are ceRNAs.
#' @param mRExp mRNA expression data, rows are samples, columns are mRNAs.
#' @param miRTarget Putative miRNA-target interactions.
#' @param minSharedmiR The minimum number of shared miRNAs between ceRNAs and mRNAs.
#' @param padjustvaluecutoff A cutoff value of adjusted p-values.
#' @param padjustmethod Adjusted method of p-values, can select one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param poscorcutoff A cutoff value of positive correlation.
#' @param senscorcutoff A cutoff value of sensitivity correlation.
#' @param null_model Pre-computed null model. Users can build the null model using "sponge_build_null_model" function in SPONGE R package.
#' @import SPONGE
#' @importFrom corpcor pcor.shrink
#' @importFrom WGCNA corAndPvalue
#' @import doParallel
#' @import foreach
#' @importFrom stats phyper
#' @importFrom stats cor.test
#' @importFrom stats p.adjust
#' @export
#' @return Matrix object: Predicted ceRNA interactions.
#' 
#' @examples 
#' # NOT RUN
#' # library(SPONGE)
#' # data(ACC)
#' # num.cores <- 2
#' # cores <- makeCluster(num.cores)
#' # registerDoParallel(cores)
#' # pre_null_model_ACC <- sponge_build_null_model(number_of_datasets = 1e+03, 
#' # number_of_samples = nrow(ACC_miRNA_Exp))
#' # ceRNet_lncR_vs_mR_ACC <- SC(ACC_miRNA_Exp, ACC_lncRNA_Exp, ACC_mRNA_Exp, miRTarget_lncR_vs_mR, null_model = pre_null_model_ACC)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Langfelder P, Horvath S. WGCNA: an R package for weighted 
#' correlation network analysis. BMC Bioinformatics. 2008, 9:559.
#' @references List M, Dehghani Amirabad A, Kostka D, Schulz MH. 
#' Large-scale inference of competing endogenous RNA networks with 
#' sparse partial correlation. Bioinformatics. 2019;35: i596–i604.
#' @references Paci P, Colombo T, Farina L. Computational analysis identifies 
#' a sponge interaction network between long non-coding RNAs and messenger RNAs 
#' in human breast cancer. BMC Syst Biol. 2014;8: 83. doi:10.1186/1752-0509-8-83
SC <- function(miRExp, 
               ceRExp, 
	             mRExp, 
	             miRTarget, 
	             minSharedmiR = 1, 
	             padjustvaluecutoff = 0.05, 
	             padjustmethod = "none", 
	             poscorcutoff = 0, 
	             senscorcutoff = 0,
	             null_model){   
         
        miRTarget <- as.matrix(miRTarget)	
        miRceR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(ceRExp))), ]
        miRmR <- miRTarget[intersect(which(miRTarget[, 1] %in% colnames(miRExp)),
	                      which(miRTarget[, 2] %in% colnames(mRExp))), ]

	      ceRSym <- unique(miRceR[, 2])
        mRSym <- unique(miRmR[, 2])
	      miRSym <- unique(c(miRceR[, 1], miRmR[, 1]))

        ceRExp_query <- ceRExp[, which(colnames(ceRExp) %in% ceRSym)]
        mRExp_query <- mRExp[, which(colnames(mRExp) %in% mRSym)]
	
        Cor.Pvalue <- corAndPvalue(ceRExp_query, mRExp_query)
        
	      index <- which(Cor.Pvalue$cor > poscorcutoff & Cor.Pvalue$p < padjustvaluecutoff, arr.ind = TRUE)
        
	# get number of cores to run	
  # cl <- makeCluster(num.cores)
  # registerDoParallel(cl)	

        Res <- foreach(i = seq_len(nrow(index)), .packages = "corpcor") %dopar% {
	    
	              Interin1 <- miRceR[which(miRceR[, 2] %in% ceRSym[index[i, 1]]), 1]
                Interin2 <- miRmR[which(miRmR[, 2] %in% mRSym[index[i, 2]]), 1]

		            M1 <- length(Interin1)
                M2 <- length(Interin2)
                M3 <- length(intersect(Interin1, Interin2))
                M4 <- length(miRSym)
                M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

                if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
                    
                    C1 <- ceRSym[index[i, 1]]
                    C2 <- mRSym[index[i, 2]]

                    ceRExpIdx <- which(colnames(ceRExp) %in% ceRSym[index[i, 1]])
                    mRExpIdx <- which(colnames(mRExp) %in% mRSym[index[i, 2]])
                    miRExpIdx <- which(colnames(miRExp) %in% intersect(Interin1, Interin2))

		                M6 <- Cor.Pvalue$cor[index[i, 1], index[i, 2]]
		                M7 <- Cor.Pvalue$p[index[i, 1], index[i, 2]]
		                M8 <- corpcor::pcor.shrink(cbind(ceRExp[, ceRExpIdx], mRExp[, mRExpIdx],
                                                     miRExp[, miRExpIdx]), verbose = FALSE)[1, 2]
                    M9 <- M6 - M8
                } else {

                C1 <- NA; C2 <- NA; M6 <- NA; M7 <- NA; M8 <- NA; M9 <- NA 
	       
                }
	       
	              tmp <- c(C1, C2, M3, M5, M6, M7, M8, M9)    
                return(tmp)
	}
	# shut down the workers
  # stopCluster(cl)
  # stopImplicitCluster()
  # unregister_dopar()

	      Res <- do.call(rbind, Res)
        
	      Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
	      Res[, 6] <- p.adjust(as.numeric(Res[, 6]), method = padjustmethod)

        Res <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 6]) < padjustvaluecutoff & as.numeric(Res[, 8]) > senscorcutoff), ]        
                
        sponge_result <- data.frame(geneA = Res[, 1], 
                                    geneB = Res[, 2], 
				                            df = as.numeric(Res[, 3]), 
				                            cor = as.numeric(Res[, 5]), 
				                            pcor = as.numeric(Res[, 7]), 
				                            mscor = as.numeric(Res[, 8]))

        sponge_result_null_model <- sponge_compute_p_values(sponge_result = sponge_result, null_model = null_model)
        
	    if(padjustmethod == "none"){
          sponge_result_fdr <- sponge_result_null_model[which(sponge_result_null_model$p.val < padjustvaluecutoff), ]
          Res_fdr <- Res[which(sponge_result_null_model$p.val < padjustvaluecutoff), ]
	    } else {
	        sponge_result_fdr <- sponge_result_null_model[which(sponge_result_null_model$p.adj < padjustvaluecutoff), ]
          Res_fdr <- Res[which(sponge_result_null_model$p.adj < padjustvaluecutoff), ]
	    }
        
	      Res_final <- cbind(Res_fdr, sponge_result_fdr[, 8])

        colnames(Res_final) <- c("sponge_1", 
	                               "sponge_2", 
				                         "Number of shared miRNAs",
				                         "p.value of shared miRNAs",
			                           "Correlation", 
				                         "p.value of correlation",
			                           "Partial correlation",
			                           "Sensitivity correlation",
				                         "p.value of sensitivity correlation")
			   
        rownames(Res_final) <- seq_len(nrow(Res_final))

return(Res_final)

}

#' Identifying ncRNA synergistic competition using Synergistic COMpetition (SCOM) method and null model.
#' 
#' @title SCOM
#' @param ceRnet Putative ceRNA interactions.
#' @param ncRExpData ncRNA expression data, rows are samples, columns are ncRNAs.
#' @param mRExpData mRNA expression data, rows are samples, columns are mRNAs.
#' @param minSharedmR The minimum number of shared competed mRNAs between ceRNAs.
#' @param padjustvaluecutoff A cutoff value of adjusted p-values.
#' @param padjustmethod Adjusted method of p-values, can select one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param poscorcutoff A cutoff value of positive correlation.
#' @param senscorcutoff A cutoff value of sensitivity correlation.
#' @param null_model Pre-computed null model. Users can build the null model using "sponge_build_null_model" function in SPONGE R package.
#' @import SPONGE
#' @importFrom corpcor pcor.shrink
#' @importFrom WGCNA corAndPvalue
#' @import doParallel
#' @import foreach
#' @export
#' @return Matrix object: Predicted ncRNA synergistic competition relationships.
#' 
#' @examples 
#' # NOT RUN
#' # library(SPONGE)
#' # data(ACC)
#' # num.cores <- 2
#' # cores <- makeCluster(num.cores)
#' # registerDoParallel(cores)
#' # pre_null_model_ACC <- sponge_build_null_model(number_of_datasets = 1e+03, 
#' # number_of_samples = nrow(ACC_miRNA_Exp))
#' # ceRNet_lncR_vs_mR_ACC <- SC(ACC_miRNA_Exp, ACC_lncRNA_Exp, ACC_mRNA_Exp, miRTarget_lncR_vs_mR, null_model = pre_null_model_ACC)
#' # ceRNet_pseudo_vs_mR_ACC <- SC(ACC_miRNA_Exp, ACC_pseudogene_Exp, ACC_mRNA_Exp, miRTarget_pseudo_vs_mR, null_model = pre_null_model_ACC)
#' # ceRNet_ACC <- rbind(ceRNet_lncR_vs_mR_ACC, ceRNet_pseudo_vs_mR_ACC)
#' # SCOMNet_ACC <- SCOM(ceRNet_ACC, cbind(ACC_lncRNA_Exp, ACC_pseudogene_Exp), ACC_mRNA_Exp, null_model = pre_null_model_ACC)
#'
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Langfelder P, Horvath S. WGCNA: an R package for weighted 
#' correlation network analysis. BMC Bioinformatics. 2008, 9:559.
#' @references List M, Dehghani Amirabad A, Kostka D, Schulz MH. 
#' Large-scale inference of competing endogenous RNA networks with 
#' sparse partial correlation. Bioinformatics. 2019;35: i596–i604.
#' @references Paci P, Colombo T, Farina L. Computational analysis identifies 
#' a sponge interaction network between long non-coding RNAs and messenger RNAs 
#' in human breast cancer. BMC Syst Biol. 2014;8: 83. doi:10.1186/1752-0509-8-83
SCOM <- function(ceRnet, 
                 ncRExpData, 
		             mRExpData, 
		             minSharedmR = 1, 
		             padjustvaluecutoff = 0.05, 
		             padjustmethod = "none", 
		             poscorcutoff = 0, 
		             senscorcutoff = 0,
		             null_model) {

        m1 <- nrow(ceRnet)
        n1 <- ncol(ceRnet)    
    
        ncRExpDataNames <- as.matrix(colnames(ncRExpData))
        mRExpDataNames <- as.matrix(colnames(mRExpData))

        ncR <- ceRnet[, 1]
        mR <- ceRnet[, 2]

        ncRSym <- unique(ncR)
        mRSym <- unique(mR)

	      ncRExpData_query <- ncRExpData[, which(colnames(ncRExpData) %in% ncRSym)]
	      mRExpData_query <- mRExpData[, which(colnames(mRExpData) %in% mRSym)]

        Cor.Pvalue <- corAndPvalue(ncRExpData_query)
        Cor.Pvalue$cor[lower.tri(Cor.Pvalue$cor, diag = TRUE)] <- 0

	      index <- which(Cor.Pvalue$cor > poscorcutoff & Cor.Pvalue$p < padjustvaluecutoff, arr.ind = TRUE)
        
        # get number of cores to run	
        # cl <- makeCluster(num.cores)
        # registerDoParallel(cl)

        Res <- foreach(i = seq_len(nrow(index)), .packages = "corpcor") %dopar% {
	    
	              Interin1 <- ceRnet[which(ceRnet[, 1] %in% ncRSym[index[i, 1]]), 2]
                Interin2 <- ceRnet[which(ceRnet[, 1] %in% ncRSym[index[i, 2]]), 2]
                
                M1 <- length(Interin1)
                M2 <- length(Interin2)
                M3 <- length(intersect(Interin1, Interin2))
                M4 <- length(mRSym)
                M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

                if (M3 >= minSharedmR & M5 < padjustvaluecutoff) {

                    C1 <- ncRSym[index[i, 1]]
                    C2 <- ncRSym[index[i, 2]]

                    ncRExpIdx1 <- which(ncRExpDataNames %in% ncRSym[index[i, 1]])
                    ncRExpIdx2 <- which(ncRExpDataNames %in% ncRSym[index[i, 2]])
                    mRExpIdx <- which(mRExpDataNames %in% intersect(Interin1, Interin2))
                    
                    # Calculate sensitivity correlation of each ncRNA-ncRNA pair
                    M6 <- cor.test(ncRExpData[, ncRExpIdx1], ncRExpData[, ncRExpIdx2])$estimate
                    M7 <- cor.test(ncRExpData[, ncRExpIdx1], ncRExpData[, ncRExpIdx2])$p.value
		                M8 <- corpcor::pcor.shrink(cbind(ncRExpData[, ncRExpIdx1], ncRExpData[, ncRExpIdx2],
                                                     mRExpData[, mRExpIdx]), verbose = FALSE)[1, 2]
		                M9 <- M6 - M8
                } else {

                C1 <- NA; C2 <- NA; M6 <- NA; M7 <- NA; M8 <- NA; M9 <- NA 
	       
                }
	       
	              tmp <- c(C1, C2, M3, M5, M6, M7, M8, M9)    
                return(tmp)
        }        
        # shut down the workers
        # stopCluster(cl)
        # stopImplicitCluster()
	      # unregister_dopar()
        
	      Res <- do.call(rbind, Res)
        
        # Extract ncRNA-ncRNA pairs with sensitivity correlation more than senscorcutoff.
        Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
	      Res[, 6] <- p.adjust(as.numeric(Res[, 6]), method = padjustmethod)

        Res <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 6]) < padjustvaluecutoff & as.numeric(Res[, 8]) > senscorcutoff), ]         
        
        sponge_result <- data.frame(geneA = Res[, 1], 
                                    geneB = Res[, 2], 
				                            df = as.numeric(Res[, 3]), 
				                            cor = as.numeric(Res[, 5]), 
				                            pcor = as.numeric(Res[, 7]), 
				                            mscor = as.numeric(Res[, 8]))

        sponge_result_null_model <- sponge_compute_p_values(sponge_result = sponge_result, null_model = null_model)
        
	    if(padjustmethod == "none"){
          sponge_result_fdr <- sponge_result_null_model[which(sponge_result_null_model$p.val < padjustvaluecutoff), ]
          Res_fdr <- Res[which(sponge_result_null_model$p.val < padjustvaluecutoff), ]
	        Res_final <- cbind(Res_fdr, sponge_result_fdr[, 7])
	    } else {
	        sponge_result_fdr <- sponge_result_null_model[which(sponge_result_null_model$p.adj < padjustvaluecutoff), ]
          Res_fdr <- Res[which(sponge_result_null_model$p.adj < padjustvaluecutoff), ]
	        Res_final <- cbind(Res_fdr, sponge_result_fdr[, 8])
	    }	

          colnames(Res_final) <- c("ncRNA_1", 
	                                 "ncRNA_2", 
		                               "Number of shared mRNAs",
			                             "p.value of shared mRNAs",
			                             "Correlation", 
		                               "p.value of correlation",
			                             "Partial correlation",
			                             "Sensitivity correlation",
		                               "p.value of sensitivity correlation")
			   
          rownames(Res_final) <- seq_len(nrow(Res_final))

return(Res_final)

}

#' Calculating the significance p-value of topological characteristics in biological networks.
#'
#' @title Random_net_parallel
#' @param obser_path Observed characteristics path length of a biological network.
#' @param obser_density Observed density of a biological network.
#' @param nodes.num The number of nodes in a biological network.
#' @param edges.num The number of edges in a biological network.
#' @param perm The number of permutations for generating random networks.
#' @param directed A logical value, FALSE for undirected network and TRUE for directed network.
#' @param num.cores The number of cores used for parallel calculation.
#' @import igraph
#' @import doParallel
#' @import foreach
#' @import parallel
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @export
#' @return A vector: The significance -log10(p-value) of characteristics path length and density of a biological network. 
#' 
#' @examples 
#' # NOT RUN
#' # library(SPONGE)
#' # data(ACC)
#' # num.cores <- 2
#' # cores <- makeCluster(num.cores)
#' # registerDoParallel(cores)
#' # pre_null_model_ACC <- sponge_build_null_model(number_of_datasets = 1e+03, 
#' # number_of_samples = nrow(ACC_miRNA_Exp))
#' # ceRNet_lncR_vs_mR_ACC <- SC(ACC_miRNA_Exp, ACC_lncRNA_Exp, ACC_mRNA_Exp, miRTarget_lncR_vs_mR, null_model = pre_null_model_ACC)
#' # ceRNet_lncR_vs_mR_ACC_graph <- make_graph(c(t(ceRNet_lncR_vs_mR_ACC[, 1:2])), directed = FALSE)
#' # ceRNet_lncR_vs_mR_ACC_path <- mean_distance(ceRNet_lncR_vs_mR_ACC_graph)
#' # ceRNet_lncR_vs_mR_ACC_density <- edge_density(ceRNet_lncR_vs_mR_ACC_graph)
#' # ceRNet_lncR_vs_mR_ACC_random <- Random_net_parallel(ceRNet_lncR_vs_mR_ACC_path, ceRNet_lncR_vs_mR_ACC_density, length(V(ceRNet_lncR_vs_mR_ACC_graph)), length(E(ceRNet_lncR_vs_mR_ACC_graph)))
#' # ceRNet_lncR_vs_mR_ACC_path_p <- ceRNet_lncR_vs_mR_ACC_random[1]
#' # ceRNet_lncR_vs_mR_ACC_density_p <- ceRNet_lncR_vs_mR_ACC_random[2]
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Csardi G, Nepusz T. The igraph software package for complex network research. 
#' InterJournal, Complex Systems. 2006, 1695.
Random_net_parallel <- function(obser_path, 
                                obser_density, 
                                nodes.num, 
                                edges.num, 
                                perm = 100, 
                                directed = FALSE, 
                                num.cores = 32) {
    
    set.seed(123)
    cores <- makeCluster(num.cores)
    registerDoParallel(cores)

    res <- foreach(i = seq(perm), .packages = "igraph") %dopar% {    
        g <- sample_pa(n = nodes.num, m = edges.num, directed = directed)
        g <- delete_edges(g, sample(1:gsize(g), size = gsize(g) - edges.num))    
        tmp_path <- mean_distance(g)    
        tmp_density <- edge_density(g)    
    
        return(c(tmp_path, tmp_density))
    }
    
    res <- do.call(rbind, res)    
    pvalue_path <- -log10(pnorm(obser_path, mean(res[, 1]), sd(res[, 1]), lower.tail = TRUE) + eps(1))    
    pvalue_density <- -log10(pnorm(obser_density, mean(res[, 2]), sd(res[, 2]), lower.tail = FALSE) + eps(1))    

    return(c(pvalue_path, pvalue_density))
}

#' Calculating similarity matrix between two list of networks
#' 
#' @title Sim.ceRNet
#' @param net1 List object, the first list of network.
#' @param net2 List object, the second list of network.
#' @param directed Logical value, network directed (TRUE) or undirected (FALSE).
#' @import igraph
#' @export
#' @return Matrix object: A similarity matrix between two list of networks.   
#' @examples
#' library(igraph) 
#' net1 <- list(as_data_frame(sample_k_regular(10, 2)), as_data_frame(sample_k_regular(10, 2))) 
#' net2 <- list(as_data_frame(sample_k_regular(10, 3)), as_data_frame(sample_k_regular(10, 3))) 
#' Sim.net <- Sim.ceRNet(net1, net2)
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references  Tucker CM, Cadotte MW, Carvalho SB, et al. A guide to phylogenetic metrics for conservation, community ecology and macroecology. Biol Rev Camb Philos Soc. 2017;92(2):698-715. 
Sim.ceRNet <- function(net1, 
                       net2, 
                       directed = FALSE){

    if(class(net1)!="list" | class(net2)!="list") {
        stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){
	          net1_graph_interin <- make_graph(c(t(net1[[i]][, 1:2])), directed = directed)
            net2_graph_interin <- make_graph(c(t(net2[[j]][, 1:2])), directed = directed)
	          overlap_interin <- nrow(as_data_frame(net1_graph_interin %s% net2_graph_interin))
	          Sim[i, j] <- overlap_interin/min(nrow(net1[[i]]), nrow(net2[[j]]))
	}
    }

    return(Sim)
}

#' Calculating similarity matrix between two list of hubs.
#' 
#' @title Sim.hub
#' @param hub1 List object, the first list of hub.
#' @param hub2 List object, the second list of hub.
#' @export
#' @return Matrix object: A similarity matrix between two list of hubs. 
#' 
#' @examples
#' library(igraph) 
#' hub1 <- list(c("ncRNA1", "ncRNA2", "ncRNA3", "ncRNA4"), c("ncRNA1", "ncRNA2", "ncRNA4", "ncRNA5")) 
#' hub2 <- list(c("ncRNA1", "ncRNA2", "ncRNA5", "ncRNA6"), c("ncRNA1", "ncRNA3", "ncRNA4", "ncRNA6"))
#' Sim_hub <- Sim.hub(hub1, hub2) 
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references  Tucker CM, Cadotte MW, Carvalho SB, et al. A guide to phylogenetic metrics for conservation, community ecology and macroecology. Biol Rev Camb Philos Soc. 2017;92(2):698-715. 
Sim.hub <- function(hub1, hub2){

    if(class(hub1)!="list" | class(hub2)!="list") {
        stop("Please check your input hub! The input hub should be list object! \n")
    }

    m <- length(hub1)
    n <- length(hub2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){	    
	        overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
	        Sim[i, j] <- overlap_interin/min(length(hub1[[i]]), length(hub2[[j]]))
	}
    }

    return(Sim)
}


#' Identifying the overlap between multiple networks.
#' 
#' @title Overlap.ceRNet
#' @param net List object, the list of network.
#' @param overlap.num The minimum number of interactions existing in multiple networks.
#' @param type The overlapped interactions in overlap.num networks ("equal") or at least overlap.num networks ("least").
#' @importFrom stringr str_split_fixed
#' @export
#' @return Matrix object: The overlapped interactions.
#' 
#' @examples
#' library(igraph) 
#' net <- list(as_data_frame(sample_k_regular(10, 2)), as_data_frame(sample_k_regular(10, 3)), as_data_frame(sample_k_regular(10, 2)), as_data_frame(sample_k_regular(10, 3))) 
#' ceRNet_2 <- Overlap.ceRNet(net, overlap.num = 2, type = "least") 
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2}) 
Overlap.ceRNet <- function(net, 
                           overlap.num = 1, 
                           type = c("equal", "least")){

    if(class(net)!="list") {
        stop("Please check your input network! The input network should be list object! \n")
    }

    net.transform <- unlist(lapply(seq(net), function(i) paste(net[[i]][, 1], net[[i]][, 2], sep = " & ")))
    net.table <- table(net.transform)
    
    if(type == "least"){
        net.overlapped <- which(net.table >= overlap.num)
    } else if(type == "equal"){
        net.overlapped <- which(net.table == overlap.num)
    }
    
    overlapped <- names(net.overlapped)
    overlapped <- str_split_fixed(overlapped, " & ", 2)
    return(overlapped)
}

#' Identifying the overlap between multiple lists of hubs.
#' 
#' @title Overlap.hub
#' @param hub List object, the list of hub
#' @param overlap.num The minimum number of hubs existing in multiple lists of hubs
#' @param type The overlapped hubs in overlap.num hub lists ("equal") or at least overlap.num hub lists ("least").
#' @export
#' @return A vector: The overlapped hubs.
#' 
#' @examples 
#' hub <- list(c("ncRNA1", "ncRNA2", "ncRNA3", "ncRNA4"), c("ncRNA1", "ncRNA2", "ncRNA4", "ncRNA5"), c("ncRNA1", "ncRNA2", "ncRNA5", "ncRNA6"), c("ncRNA1", "ncRNA3", "ncRNA4", "ncRNA6")) 
#' hub_2 <- Overlap.hub(hub, overlap.num = 2, type = "least") 
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2}) 
Overlap.hub <- function(hub, 
                        overlap.num = 1, 
                        type = c("equal", "least")){

    if(class(hub)!="list") {
        stop("Please check your input hub! The input hub should be list object! \n")
    }

    hub.transform <- unlist(hub)
    hub.table <- table(hub.transform)
    
    if(type == "least"){
        hub.overlapped <- which(hub.table >= overlap.num)
    } else if(type == "equal"){
        hub.overlapped <- which(hub.table == overlap.num)
    }
    
    overlapped <- names(hub.overlapped)
    return(overlapped)
}

#' Evaluating the performance of gene list in classifying tumor types.
#' 
#' @title module.classify
#' @param Exp Gene expression data, rows are samples, columns are genes.
#' @param tumor_type Molecular subtype information of a tumor.
#' @param Modulelist A list of gene modules.
#' @param method The multi-label classification method. It also accepts the name of the method as a string.
#' @param base.algorith The base algorithm of the multi-label classification method.
#' @param cv.folds Number of folds (Default: 10).
#' @param cv.sampling The method to split the data. (Default: "stratified")
#' @param cv.seed An optional integer used to set the seed.
#' @import mldr
#' @import utiml
#' @import e1071
#' @export
#' @return Matrix object: Multi-label classification results of a tumor.
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2}
#' @references Tsoumakas G, Katakis I, Vlahavas I. Mining multi-label data. 
#' Data mining and knowledge discovery handbook Springer, Boston, MA. Springer, Boston, MA; 2009. pp. 667–685.
#' @references Rivolli A, Carvalho ACPLF de. The utiml package: multi-label 
#' classification in R. The R Journal. 2018;10: 24–37. doi:10.32614/RJ-2018-041
#' @references Chang CC, Lin CJ. LIBSVM: a library for support vector machines. 
#' ACM Transactions on Intelligent Systems and Technology. 2011;2: 1–27.
#' @references Meyer D, Dimitriadou E, Hornik K, Weingessel A, Leisch F, Chang C-C, et al. 
#' e1071: misc functions of the department of statistics, probability theory group (Formerly: E1071), 
#' TU Wien. Available: https://CRAN.R-project.org/package=e1071
module.classify <- function(Exp, 
                            tumor_type, 
                            Modulelist, 
                            method = "br", 
                            base.algorith = "SVM", 
                            cv.folds = 10, 
	                          cv.sampling = "stratified", 
                            cv.seed = 12345) {

    module_Exp <- lapply(seq_along(Modulelist), function(i) Exp[, which(colnames(Exp) %in% Modulelist[[i]])])
    unique_type <- unique(tumor_type[, 2])
    class_infor <- do.call(cbind, lapply(seq(unique_type), function(i) as.numeric(tumor_type[, 2] == unique_type[i])))
    module_classify <- list()

    for (i in seq_along(Modulelist)){        
	
        temp <- as.data.frame(cbind(module_Exp[[i]], class_infor))
	      Indices <- ncol(temp)
	      temp_mldr <- mldr_from_dataframe(temp, labelIndices = (Indices-length(unique_type)+1):Indices, name = "TEMPMLDR")
        temp_res <- cv(temp_mldr, method = method, base.algorith = base.algorith, cv.folds = cv.folds, 
	                     cv.sampling = cv.sampling, cv.seed = cv.seed)
        module_classify[[i]] <- temp_res

    }

    return(module_classify)
}

#' Evaluating the performance of ncRNA synergistic competition for classifying BRCA subtypes.
#' 
#' @title BRCA.classify
#' @param Exp The gene expression data, rows are samples, columns are genes.
#' @param BRCA_subtype Molecular subtype information of BRCA.
#' @param genename A list of gene names.
#' @param method The multi-label classification method. It also accepts the name of the method as a string.
#' @param base.algorith The base algorithm of the multi-label classification method.
#' @param cv.folds Number of folds (Default: 10).
#' @param cv.sampling The method to split the data. (Default: "stratified")
#' @param cv.seed An optional integer used to set the seed.
#' @import mldr
#' @import utiml
#' @import e1071
#' @export
#' @return Matrix object: Multi-label classification results of BRCA.
#' 
#' @author Junpeng Zhang (\url{https://www.researchgate.net/profile/Junpeng-Zhang-2})
#' @references Tsoumakas G, Katakis I, Vlahavas I. Mining multi-label data. 
#' Data mining and knowledge discovery handbook Springer, Boston, MA. Springer, Boston, MA; 2009. pp. 667–685.
#' @references Rivolli A, Carvalho ACPLF de. The utiml package: multi-label 
#' classification in R. The R Journal. 2018;10: 24–37. doi:10.32614/RJ-2018-041
#' @references Chang CC, Lin CJ. LIBSVM: a library for support vector machines. 
#' ACM Transactions on Intelligent Systems and Technology. 2011;2: 1–27.
#' @references Meyer D, Dimitriadou E, Hornik K, Weingessel A, Leisch F, Chang C-C, et al. 
#' e1071: misc functions of the department of statistics, probability theory group (Formerly: E1071), 
#' TU Wien. Available: https://CRAN.R-project.org/package=e1071
BRCA.classify <- function(Exp, 
                          BRCA_subtype, 
                          genename, 
                          method = "br", 
                          base.algorith = "SVM", 
                          cv.folds = 10, 
	                        cv.sampling = "stratified", 
                          cv.seed = 12345) {

    gene_Exp <- Exp[, which(colnames(Exp) %in% genename)]
    Basal <- as.numeric(BRCA_subtype[, 2] == "Basal")
    Her2 <- as.numeric(BRCA_subtype[, 2] == "Her2")
    LumA <- as.numeric(BRCA_subtype[, 2] == "LumA")
    LumB <- as.numeric(BRCA_subtype[, 2] == "LumB")
    Normal <- as.numeric(BRCA_subtype[, 2] == "Normal")    
    
    temp <- as.data.frame(cbind(gene_Exp, Basal, Her2, LumA, LumB, Normal))
    Indices <- ncol(temp)
    temp_mldr <- mldr_from_dataframe(temp, labelIndices = c(Indices-4, Indices-3, Indices-2, Indices-1, Indices), name = "TEMPMLDR")
    BRCA_classify <- cv(temp_mldr, method = method, base.algorith = base.algorith, cv.folds = cv.folds, 
	                      cv.sampling = cv.sampling, cv.seed = cv.seed)
    
    return(BRCA_classify)
}
