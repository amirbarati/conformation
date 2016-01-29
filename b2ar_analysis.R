library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(classInt)
library(maptools)
library(mixtools)
library(glmnet)
library(randomForest)
library(tree)
library(R.matlab)
library(rpart)
library(foreach)
library(doMC)
registerDoMC(16)
library(randomForest)

if(file.exists("/Users/evan/vsp/b2ar_analysis/conformation/analysis.R")) {
  analysis.functions <- "/Users/evan/vsp/b2ar_analysis/conformation/analysis.R"
} else {
  analysis.functions <- "/home/enf/b2ar_analysis/conformation/analysis.R"
}
source(analysis.functions)

do.docking.analysis <- function(docking, reference.docking, tica.coords.averages, pnas.coords.averages, refcoords, analysis.dir, inverse.agonists) {
  reference.docking <- data.frame(read.csv("/home/enf/b2ar_analysis/reference_docking/docking_SP/all_docking_combined_manual.csv", stringsAsFactors = F, row.names=1))
  
  print("starting docking analysis ")
  print(inverse.agonists)
  aggregate.docking <- as.data.frame(compute.aggregate.docking(docking, inverse.agonists))
  docking.averages <- cluster.averages.new(aggregate.docking)
  write.table(docking.averages, paste(analysis.dir, "docking.averages.csv", sep="/"))
  
  docking.tica <- combine.dfs(docking.averages, tica.coords.averages)
  docking.tica <- docking.tica[order(docking.tica[,1]*-1.0),]
  colnames(docking.tica)[1] <- "rel.docking.score"
  
  rel.docking.tica <- docking.tica
  rel.docking.tica[,1] <- convert.docking.to.rel.probs(rel.docking.tica[,1])
  rel.docking.tica <- as.data.frame(rel.docking.tica)
  print(rel.docking.tica[1:2,])
  write.table(rel.docking.tica, paste(analysis.dir, "rel.docking.tica.csv", sep = "/"))
  
  #x <- model.matrix()
  rel.docking.tica.logistic <- glm(rel.docking.score~.,data=rel.docking.tica,family=binomial)
  logistic.cv <- cv.glm(data=rel.docking.tica,glmfit = rel.docking.tica.logistic,K=5)
  write.table(x = rel.docking.tica.logistic$coef, file = paste(analysis.dir, "/rel_docking_tica_logistic_coef.csv",sep=""))
  write.table(x = logistic.cv$delta, file = paste(analysis.dir, "/rel_docking_tica_logistic_cv_delta.csv",sep=""))
  
  tica.transformed <- tica.coords.averages
  intercept = summary(rel.docking.tica.logistic)$coefficients[1,1]
  coefficients = summary(rel.docking.tica.logistic)$coefficients[-1,1]
  tica.transformed <- logistic.transform(tica.transformed,intercept,coefficients)
  tica.transformed.predict <- predict.glm(rel.docking.tica.logistic, newdata=rel.docking.tica,type="response")
  #tica.unclustered.transformed <- logistic.transform(tica.coords,intercept,coefficients)
  #plot(hexbin(tica.transformed,tica.coords.averages[,1]))
  #plot(hexbin(tica.coords[,6], tica.unclustered.transformed))
  #plot.colmap(pnas.coords.averages[,3,drop=F],as.data.frame(cbind(tica.coords.averages[,1], tica.transformed)), refcoords[,coords], "tIC1_vs_transformed_tIC_w_NPxxY_color", analysis.dir)
  #plot.colmap(as.data.frame(tica.transformed),pnas.coords.averages[,c(1,3)], refcoords[,coords], "A priori Reaction Coordinates vs. transformed tIC 1", analysis.dir)
  plot(tica.transformed,pnas.coords.averages[,3])
  plot(tica.coords.averages[,"tIC.9"],pnas.coords.averages[,3])
  
  #plot.docking.vs.reaction.coord(docking.averages, pnas.coords.averages, analysis.dir)
  
  active.rows <- find_active(pnas.coords.averages)
  test_method(docking.averages, active.rows, "Aggregate Docking Score, Five Agonists, Two Inverse Agonists", analysis.dir)
  plot.docking.vs.reference(docking.averages, reference.docking, analysis.dir)
  
  #docking.subset <- docking
  #docking.subset <- docking[,which(colnames(docking) %in% inverse.agonists)]
  #docking.subset <- docking[,c("X3p0g_lig"),drop=F]
  

  coords <- c("tm6_tm3_dist", "npxxy_rmsd_active")
  #plot.colmap(docking.averages, pnas.coords.averages[,coords, drop=F], refcoords[,coords], "Docking Score All Clusters vs TM6_TM3 dist and NPxxY RMSD to Inactive", analysis.dir, top=100)
  active.docking <- combine.dfs(docking.averages, as.data.frame(active.rows))
  active.docking.sorted <- active.docking[order(-1.0*active.docking[,1]),]
  
  #Lasso with setting all probabilities greater than 0.5 to 1 and all probabilities less than 0.5 to 0.0
  X <- model.matrix(rel.docking.score ~ . -1, data = rel.docking.tica)
  y <- rel.docking.tica$rel.docking.score
  low.rows <- (y < (mean(y)+0*sd(y)))
  y[low.rows] <- 0.0
  y[!low.rows] <- 1.0
  fit.lasso <- glmnet(x = X, y = y,family="binomial")
  
  pdf(paste(analysis.dir, "/lasso_logistic_cutoff_0pt5.pdf", sep= ""))
  plot(fit.lasso, xvar = "lambda", label=TRUE)
  dev.off()
  
  cv.lasso = cv.glmnet(X, y,nfolds=5)
  min(cv.lasso$cvm)
  pdf(paste(analysis.dir, "/lasso_logistic_cutoff_0pt5_cv_vs_lambda.pdf", sep= ""))
  plot(cv.lasso)
  dev.off()
  
  #Logistic by doing 1. logit transform on rel docking dist and then 2. fitting linear function
  rel.docking.logit <- rel.docking.tica 
  rel.docking.logit[,1] <- log(rel.docking.logit[,1] / (1 - rel.docking.logit[,1]))
  
  X <- model.matrix(rel.docking.score ~ . -1, data = rel.docking.logit)
  y <- rel.docking.tica$rel.docking.score
  
  fit.logit.lasso <- glmnet(x = X, y = y)
  
  pdf(paste(analysis.dir, "/lasso_logistic_logit.pdf", sep= ""))
  plot(fit.logit.lasso, xvar = "lambda", label=TRUE)
  dev.off()
  
  cv.lasso = cv.glmnet(X, y,nfolds=5)
  min(cv.lasso$cvm)
  pdf(paste(analysis.dir, "/lasso_logistic_cutoff_0pt5_cv_vs_lambda.pdf", sep= ""))
  plot(cv.lasso)
  dev.off()
  
  
  #summary(rel.docking.tica.logistic)
  summary(fit.lasso)
  
  #backward subset selection:
  #set.seed(11)
  #nfolds <- 5
  #folds <- sample(rep(1:nfolds),length=nrow(rel.docking.tica))
  #cv.errors <- matrix(NA, nfolds, ncol(tica.coords.averages))
  
  #for (k in 1:5) {
  #  b
  #}
  #Xy <- as.data.frame(cbind(x,y))
  #bestBIC <- bestglm(Xy, IC="BIC",family= quasibinomial)
  
  
}

do.analysis <- function(tica, analysis.dir, pnas.coords.csv, tica.coords.csv, features.dir = "", docking.csv="", sasa.csv="", docking.aggregated.csv="") {
 #print("Conducting Analysis in R")
 #print(tica)
 #print(analysis.dir)
 #print(pnas.coords.csv)
 #print(tica.coords.csv)
 #print(features.dir)
  inverse.agonists <- c("s.atenolol", "s.carazolol")

  #base <- "/Users/Evan/vsp/b2ar_analysis"
  #base <- "/Users/Evan/vsp/b2ar_analysis/exacycle_data"
  #tica <- "/tICA_t5_n_components10all_residues_under_cutoff1nm_regularizationdefault"
  #tica <- "/tICA_t5_n_components10_skip5_switches_pp_npxx_contact_cutoff10000nm_regularization2pt0"
  #tica <- "/tICA_t5_n_components10skip5_switches_pp_npxx_ser_cutoff10000nm_regularization_0pt5/"
  #tica <- "/tICA_t5_n_components10skip5_switches_pp_npxx_ser_cutoff10000nm_regularizationdefault"
  #tica <- "/tICA_t5_n_components10_skip5_switches_pp_npxx_contact"
  #tica <- "/tICA_t5_n_components10_skip5_switches_pp_npxx_contact/ktICA_n_components5_random_specified"
  #tica <- "/tICA_t10_n_components10_skip5_switches_pp_npxx_contact"
  #tica <- "/tICA_t10_n_components5_switches_npxx_tm6_bp"
  #tica <- "tICA_t10_n_components"
  #tica <- paste(base, tica, sep = "")
  #clusters <- 1000
  #exacycle <- ""
 # method <- "random"
  #analysis.dir <- paste(tica, "/analysis_n_clusters", clusters, "_", method, sep = "")

  #pnas.coords.csv <- paste(analysis.dir, "/pnas_coords_new.csv", sep = "")
  #tica.coords.csv <- paste(analysis.dir, "/tica_coords.csv", sep = "")
  
  #docking.csv <- paste(analysis.dir, "/all_docking_combined.csv", sep = "")
  #docking.aggregated.csv <- paste(analysis.dir, "/aggregate_docking.csv", sep = "")
  #sasa.csv <- paste(analysis.dir, "/sasa_bp.csv", sep = "")
  #pnas.coords.all_csv <- "/Users/evan/scratch_enf/b2ar_analysis/exacycle_data/all_pnas_features/pnas_all_coords.csv"
  #pnas.coords.all <- data.frame(read.csv(pnas.coords.all_csv, stringsAsFactors = F, row.names=1))
  #badrows <- which(pnas.coords.all[,2] > 4.0 | pnas.coords.all[,3] > 4.0 | pnas.coords.all[,4] > 4.0 | pnas.coords.all[,5] > 4.0)
  #if(length(badrows) > 0) {
  #  pnas.coords.all <- pnas.coords.all[-badrows,]
  #}
  #colnames(pnas.coords.all) <- colnames(pnas.coords)
  #pnas.coords.all["tm6_tm3_dist"] <- 7.14 * pnas.coords.all["tm6_tm3_dist"]
  
  
  pnas.coords <- data.frame(read.csv(pnas.coords.csv, stringsAsFactors = F, row.names=1))[,c(1,2,3,4,5)]
  pnas.coords["tm6_tm3_dist"] <- 7.14 * pnas.coords["tm6_tm3_dist"]
  tica.coords <- data.frame(read.csv(tica.coords.csv, stringsAsFactors = F, row.names=1))
  for (i in 1:dim(tica.coords)[2]) {
    colnames(tica.coords)[i] <- paste("tIC.", i, sep="")
  }
  badrows <- which(pnas.coords[,2] > 4.0 | pnas.coords[,3] > 4.0 | pnas.coords[,4] > 4.0 | pnas.coords[,5] > 4.0)
  if(length(badrows) > 0) {
    pnas.coords <- pnas.coords[-badrows,]
    tica.coords <- tica.coords[-badrows,]
  }
  dim(pnas.coords)
  pnas.coords[1:10,]
  
  pnas.coords.averages <- cluster.averages.new(pnas.coords)
  write.table(pnas.coords.averages, file = paste(analysis.dir, "/pnas_coords_averages.csv", sep=""))
  tica.coords.averages <- cluster.averages.new(tica.coords)
  write.table(tica.coords.averages, file = paste(analysis.dir, "/tica_coords_averages.csv", sep=""))
  
  refcoords <- read.csv(file = "/home/enf/b2ar_analysis/reference_receptors/ref_coords.csv", row.names = 1)
  colnames(refcoords) <- colnames(pnas.coords)
  coords <- c("tm6_tm3_dist", "npxxy_rmsd_active")
  if(1==2) {
  if(file.exists(paste(tica, "/reference_receptors/refcoords.csv", sep=""))) {
    reftica <- read.csv(file = paste(tica, "/reference_receptors/refcoords.csv", sep=""),header=F)
    for(i in 1:dim(tica.coords.averages)[2]) {
      plot.colmap(tica.coords.averages[,i,drop=F],pnas.coords.averages[,coords,drop=F], refcoords[,coords], paste("A priori Reaction Coordinates vs. tIC ", i,sep=""), analysis.dir)
      ggplot.colmap(tica.coords.averages[,i,drop=F],pnas.coords.averages[,coords,drop=F], refcoords[,coords], reftica[,i,drop=F], paste("A priori Reaction Coordinates vs. tIC ", i,sep=""), analysis.dir)
    }
  } else {
    for(i in 1:dim(tica.coords.averages)[2]) {
      plot.colmap(tica.coords.averages[,i,drop=F],pnas.coords.averages[,coords,drop=F], refcoords[,coords], paste("A priori Reaction Coordinates vs. tIC ", i,sep=""), analysis.dir)
      ggplot.colmap(tica.coords.averages[,i,drop=F],pnas.coords.averages[,coords,drop=F], refcoords[,coords], "", paste("A priori Reaction Coordinates vs. tIC ", i,sep=""), analysis.dir)
    }
  }
  }
  


  
  pnas.tica <- combine.dfs(pnas.coords.averages, tica.coords.averages)[,c(1,3,5,6,7,8,9,10,11,12,13,14,15),drop=F]
  pdf(paste(analysis.dir, "/pnas_tica_pairs.pdf",sep=""))
  pairs(pnas.tica,pch='.')
  dev.off()
  
  #colnames(pnas.coords.all) <- colnames(pnas.coords)
  #pnas.coords.all["tm6_tm3_dist"] <- 7.14 * pnas.coords.all["tm6_tm3_dist"]
  
  #Find which residues are most important in each tIC

  if(length(grep("ktICA", tica)) == 100000) {
    tic.residues <- get.tic.residues.pairs(analysis.dir = analysis.dir, features.dir = features.dir, tica = tica)
    #print(tic.residues[1:10,])
    residue.names <- combine.columns(tic.residues[,dim(tic.residues)[2]-1], tic.residues[,dim(tic.residues)[2]])
    rownames(tic.residues) <- residue.names
    tic.residues <- tic.residues[,c(1:(dim(tic.residues)[2]-2))]
    tic.residue.sorted <- sort.by.column(tic.residues)
    for(i in 1:dim(tic.residues)[2]) {
      pdf(paste(analysis.dir, "/barplot_residue_pair_coefficients_tIC", i, ".pdf", sep=""),width=20,height=5)
      df <- tic.residue.sorted[[i]]
      df1 <- df[1:100,,drop=F]
      #df2 <- as.data.frame(matrix(0, nrow=20, ncol=dim(df)[2]))
      #colnames(df2) <- colnames(df1)
      #df3 <- df[(dim(df)[1]-20) : (dim(df)[1]),,drop=F]
      #df <- rbind(df1,df2,df3)
      #print(df)
      barplot(df1[,i],las=2, main = paste("tIC",i, " residue pair coefficients", sep=""))
      dev.off()
    }
    
    
    tic.residues <- get.tic.residues(analysis.dir = analysis.dir, features.dir = features.dir, tica = tica)
    #print("tic.residues:")
    #print(tic.residues)[1:10,]
    tic.residue.averages <- cluster.averages.new(tic.residues)
    tic.residue.averages.sorted <- sort.by.column(tic.residue.averages)
    for(i in 1:dim(tic.residue.averages)[2]) {
      pdf(paste(analysis.dir, "/barplot_average_residue_coefficients_tIC", i, ".pdf", sep=""),width=20,height=5)
      barplot(tic.residue.averages.sorted[[i]][,i],las=2, main = paste("tIC",i, " average residue coefficients", sep=""))
      dev.off()
    }
    
    tic.residues.95th.percentile <- cluster.percentiles.new(tic.residues, 0.95)
    tic.residues.95th.percentile.sorted <- sort.by.column(tic.residues.95th.percentile)
    for(i in 1:dim(tic.residue.averages)[2]) {
      pdf(paste(analysis.dir, "/barplot_95th_percentile_residue_coefficients_tIC", i, ".pdf", sep=""))
      barplot(tic.residues.95th.percentile.sorted[[i]][,i],las=2, main = paste("tIC",i, " 95th percentile residue coefficients", sep=""))
      dev.off()
    }
    
    tic.residues <- get.tic.residues(analysis.dir = analysis.dir, features.dir = features.dir, tica = tica)
    helix <- which.helix(tic.residues, only.helices=TRUE)
   #print("Helix:")
   #print(dim(helix))
   #print(helix[1:10,])
    helix.averages <- cluster.averages.new(helix)
    helix.sorted <- sort.by.column(helix.averages)
    for(i in 1:dim(tic.residue.averages)[2]) {
      pdf(paste(analysis.dir, "/barplot_average_helix_coefficient_tIC", i, ".pdf", sep=""))
      barplot(helix.sorted[[i]][,i],las=2, main = paste("tIC",i, " average helix coefficient", sep=""))
      dev.off()
    }
    
    tic.residues <- get.tic.residues.pairs(analysis.dir = analysis.dir, features.dir = features.dir, tica = tica)
    interhelix <- interhelix.distance(tic.residues)
    interhelix.averages <- cluster.averages.new(interhelix)
    interhelix.averages.sorted <- sort.by.column(interhelix.averages)
    for(i in 1:dim(interhelix.averages)[2]) {
      pdf(paste(analysis.dir, "/barplot_average_inter_helix_coefficient", i, ".pdf", sep=""))
      barplot(interhelix.averages.sorted[[i]][,i],las=2, main = paste("tIC",i, " average interhelix distance coefficient", sep=""))
      dev.off()
    }
    
    tic.residues <- get.tic.residues.pairs(analysis.dir = analysis.dir, features.dir = features.dir, tica = tica)
    interhelix <- interhelix.distance(tic.residues, only.helices=TRUE)
   #print("interhelix:")
   #print(dim(interhelix))
   #print(interhelix[1:10,])
    interhelix.averages <- cluster.averages.new(interhelix)
    interhelix.averages.sorted <- sort.by.column(interhelix.averages)
    for(i in 1:dim(interhelix.averages)[2]) {
      pdf(paste(analysis.dir, "/barplot_average_inter_helix_only_coefficient", i, ".pdf", sep=""))
      barplot(interhelix.averages.sorted[[i]][,i],las=2, main = paste("tIC",i, " average interhelix distance coefficient", sep=""))
      dev.off()
    }
    
    tic.residues <- get.tic.residues.pairs(analysis.dir = analysis.dir, features.dir = features.dir, tica = tica)
    interhelix <- interhelix.distance(tic.residues, only.helices=TRUE)
    interhelix.averages <- cluster.percentiles.new(interhelix,0.95)
    interhelix.averages.sorted <- sort.by.column(interhelix.averages)
    for(i in 1:dim(interhelix.averages)[2]) {
      pdf(paste(analysis.dir, "/barplot_percentile_0pt95_inter_helix_only_coefficient", i, ".pdf", sep=""))
      barplot(interhelix.averages.sorted[[i]][,i],las=2, main = paste("tIC",i, " 95th percentile interhelix distance coefficient", sep=""))
      dev.off()
    }
  }
  
  #tic.residues.sum <- cluster.sums(tic.residues)
  #tic.residues.sum.sorted <- sort.by.column(tic.residues.sum)

  active.rows <- find_active(pnas.coords.averages)
  write.table(t(as.data.frame(names(active.rows[active.rows==T]))), file = paste(analysis.dir, "/", "active_clusters.csv", sep=""), row.names = F, col.names = F, sep = ",")
  
  intermediate.rows <- find.intermediate(pnas.coords.averages)
  write.table(t(as.data.frame(names(intermediate.rows[intermediate.rows==T]))), file = paste(analysis.dir, "/", "intermediate_clusters.csv", sep=""), row.names = F, col.names = F, sep = ",")
  
  inactive.rows <- find.inactive(pnas.coords.averages)
  write.table(t(as.data.frame(names(inactive.rows[inactive.rows==T]))), file = paste(analysis.dir, "/", "inactive_clusters.csv", sep=""), row.names = F, col.names = F, sep = ",")
  

  #do.docking.analysis(docking, reference.docking, tica.coords.averages, pnas.coords.averages)
  
 print(docking.csv)
  if(file.exists(docking.csv)){
    print("Doing docking analysis")
    print(inverse.agonists)
    docking <- data.frame(read.csv(docking.csv, stringsAsFactors = F, row.names=1))
    refcoords <- read.csv(file = "/home/enf/b2ar_analysis/reference_receptors/ref_coords.csv", row.names = 1)
    colnames(refcoords) <- colnames(pnas.coords)
    do.docking.analysis(docking, reference.docking, tica.coords.averages, pnas.coords.averages, refcoords, analysis.dir, inverse.agonists)
    #docking.aggregated <- data.frame(read.csv(docking.aggregated.csv, stringsAsFactors = F, row.names=1))
    #docking.averages <- cluster.averages.new(aggregate.docking)
    #rel.docking.tica.linear <- glm(rel.docking.score~.,data=rel.docking.tica)
    #summary(rel.docking.tica.linear)
    #summary(rel.docking.tica.logistic)
   
  }
  #sasa <- data.frame(read.csv(sasa.csv, stringsAsFactors = F, row.names=1))
  

  #plot_coords_and_docking(pnas.coords.docking.sorted[1:50,], pnas.coords.docking.sorted, refcoords)

  #all_pnas_rows <- find_active(pnas.coords.all)
  #all_active_rows <- all_pnas_rows[all_pnas_rows == T]

  #

  
  #pnas.rows <- find_active(pnas.coords)
  #docking.aggregated <- as.data.frame(apply(docking.averages[,c(1,2,3,4),drop=F],1,mean))
  

  #test_method(docking.averages, pnas.rows)
  #all_pnas_rows <- find_active(pnas.coords.all)
  #all_active_rows <- all_pnas_rows[all_pnas_rows == T]
  
  #docking.csv <- "/Users/evan/biox3/b2ar_analysis_sherlock_all/b2ar_analysis/tICA_t10_n_components5_switches_npxx_tm6_bp/analysis_n_clusters1000/aggregate_docking_joined.csv"
  
  
  
  #0.779 for cutoff of 0.3
  #0.818 for cutoff of 0.25
  #0.948 for cutoff of 0.2
  
  #print(length(active_rows))
  #print(length(active_rows)/dim(pnas.coords.averages)[1])
  #print(length(all_active_rows)/(length(all_pnas_rows)))
  
  
  
  #pnas_clusters_rows <- find_active(pnas.coords.averages)
  #active_clusters_rows <- pnas_clusters_rows[pnas_clusters_rows == T]
  #inactive_clusters_rows <- pnas_clusters_rows[pnas_clusters_rows == F]
  
  #pnas_rows <- find_active(pnas.coords)
  #active_rows <- pnas_rows[pnas_rows == T]
  #inactive_rows <- pnas_rows[pnas_rows == F]
  
 #print(length(active.rows[active.rows==T]))
  write(length(active.rows[active.rows==T]), file = paste(analysis.dir, "/num_active_clusters.txt",sep=""))
  #print(length(active_rows)/dim(pnas.coords)[1])
  #print(length(all_active_rows)/(length(all_pnas_rows)))
  
  #print(length(active_clusters_rows))
  #print(length(active_clusters_rows)/dim(pnas.coords.averages)[1])
  #print(length(all_active_rows)/(length(all_pnas_rows)))
  
  #hist(pnas.coords[,1], breaks=20)
  #hist(tica.coords[,2], breaks=20)
}

#CITE: Lecture 20, Mixture Examples and Complements, 36-402, Adavnced Data Analysis, 5 APril 2011

BIC.mix <- function(out) {
    return(-2*out$loglik+log(length(out$x))*(length(unlist(out[2:4]))-1))
}

plot.mixture <- function(tic, best.model, n.comp, filename) {
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    gp <- ggplot(tic, aes_string(x = colnames(tic)[1])) + geom_histogram(aes(y= ..density..), binwidth=0.1) + geom_density()

    ggplot.normal.components <- function(component.number, mixture,...) {
        comp.fun <- function(x) {
         return(mixture$lambda[component.number] *
              dnorm(x,mean=mixture$mu[component.number],
                    sd=mixture$sigma[component.number]))
     }
        return(stat_function(comp.fun))
    }

    mixture= best.model
    MyFun <- function(x, component.number) {
        return(mixture$lambda[component.number] *
              dnorm(x,mean=mixture$mu[component.number],
                    sd=mixture$sigma[component.number]))
    }
    coeflines = list()

    for(component.number in 1:n.comp) {
        print(component.number)
        coeflines[[component.number]] = stat_function(fun=MyFun, n =1000, args = list(component.number=component.number), col = cbPalette[component.number])
    }
    gp <- gp + coeflines #+ scale_colour_manual("Lgend title", values=cbPalette[1:5])
    ggsave(filename = filename, plot = gp)
}

compute.tIC.mixture.model <- function(tIC, j, save.dir, max_components = 10, num.repeats=5) {
  min.components <- rep(NA, num.repeats)
  for(r in 1:num.repeats) {
    bics <- rep(NA, max_components)
    loglikes <- rep(NA, max_components)
    n <- length(tIC)
    data.points <- 1:n
    data.points <- sample(data.points)
    train <- data.points[1:floor(n/2)]
    test <- data.points[-(1:floor(n/2))]

    loglikes <- vector(length=max_components)
    mu0 <- mean(tIC[train])
    sd0 <- sd(tIC[train] * sqrt(floor(n/2)-1)/floor(n/2))

    mixtures = list()
    mixtures[[1]] <- list(x=tIC[train], mu=mu0,sigma=sd0,lambda=1,
                    loglik=sum(dnorm(tIC[test],mu0,sd0,log=TRUE)))
    bics[1] <- BIC.mix(mixtures[[1]])
    loglikes[1] <- sum(dnorm(tIC[test],mu0,sd0,log=TRUE))
    candidate.component.numbers <- 2:max_components
    for(k in candidate.component.numbers) {
      print(paste("fitting model for ", k, " components", sep=""))
      mixtures[[k]] <- tryCatch({
        if(k == 2) { 
          ini_mu = c(mean(tIC) - 2 * sd(tIC), mean(tIC) + 2 * sd(tIC))
        } else {
          ini_mu=seq(mean(tIC) - 2*sd(tIC), mean(tIC) + 2*sd(tIC),4*sd(tIC)/(k-1))
        }
        model <- normalmixEM(tIC[train], k=k, maxit=400,epsilon=1e-2, arbvar=FALSE, sd.constr = rep('a',k), mu=ini_mu, arbmean=TRUE)
        bics[k] <- BIC.mix(model)
        loglikes[k] <- loglike.normalmix(tIC[test],mixture=model)
        model
      }, error = function(e) {
        print("Error in fitting normal mix EM!")
        bics[k] <- 10000000.0
        loglikes[k] <- -10000000.0
      }, finally = {
      })
    }
    num.components <- min(which.min(bics),which.max(loglikes))
    min.components[r] <- num.components
  }

  print(min.components)
  k <- min(min.components)
  if(k == 2) { 
    ini_mu = c(mean(tIC) - 2 * sd(tIC), mean(tIC) + 2 * sd(tIC))
  } else {
    ini_mu=seq(mean(tIC) - 2*sd(tIC), mean(tIC) + 2*sd(tIC),4*sd(tIC)/(k-1))
  }
  final.model <- normalmixEM(tIC, k=k, maxit=400, epsilon=1e-2,arbvar=FALSE, sd.constr = rep('a',k), mu=ini_mu, arbmean=TRUE)
  final.means <- final.model$mu
  final.classes <- apply(final.model$posterior,1,which.max)

  save(final.model, file=paste(save.dir, "/tIC_", j, "_gmm_R.RData", sep=""))

  final.model.list = list()
  final.model.list$mu = final.means
  final.model.list$classes = final.classes

  return(final.model.list)

}

tIC.mixture.models <- function(tics = "", tica.csv, dir) {
  tica.csv = "/home/enf/b2ar_analysis/tICA_t5_n_components25all_tm_residues_under_cutoff1nm_regularization0pt5/analysis_n_clusters1000_random/tica_coords.csv"
  dir = "/home/enf/b2ar_analysis/tICA_t5_n_components25all_tm_residues_under_cutoff1nm_regularization0pt5/analysis_n_clusters1000_random/"
  tics = ""
  if(tics == "") {
      tics = read.csv(tica.csv, stringsAsFactors = F, header = T, row.names = 1)
  }
  
  #rows <- seq(1,dim(tics)[1],1)
  #random.rows <- sample(rows, 100000)
  #tics <- tics[random.rows,]
  best.components <- rep(NA,dim(tics)[2])
  tica.classes <- tics
  tic.names <- rep(NA, dim(tics)[2])
  for(i in 1:dim(tics)[2]) {
    tic.names[i] <- paste("tIC.", i, sep="")
  }
  colnames(tics) <- tic.names


  for(i in 1:dim(tics)[2]) {
    bics <- rep(NA,5)
      print(paste("Analyzing tIC ", i))

    #pdf(paste(dir, "/tIC.", i, "_histogram.pdf", sep=""))
    #hist(tics[,i],breaks=100, xlab = paste("tIC.", i, sep=""))
    #dev.off()
    if(1==1) {
    n <- dim(tics)[1]
    data.points <- 1:n
    data.points <- sample(data.points)
    train <- 1:n
    #train <- data.points[1:floor(n/2)]
    test <- data.points[-(1:floor(n/2))]
    candidate.component.numbers <- 2:5
    loglikes <- vector(length=1+length(candidate.component.numbers))
    mu <- mean(tics[train,i])
    sigma <- sd(tics[train,i] * sqrt(n-1)/n)
    mixtures = list()
    mixtures[[1]] <- list(x=tics[,i],mu=mean(tics[,i]),sigma=sd(tics[,i]),lambda=1,
              loglik=sum(dnorm(tics[,i],mean(tics[,i]),
                               sd(tics[,i]),log=TRUE)))
    bics[1] <- BIC.mix(mixtures[[1]])
    #loglikes[1] <- sum(dnorm(tics[test,i],mu,sigma,log=TRUE))
    #training.loglikes <- vector(length=1+length(candidate.component.numbers))
    #training.loglikes[1] <- sum(dnorm(tics[train,i], mu, sigma, log=TRUE))
    for(k in candidate.component.numbers) {
      print(paste("fitting model for ", k, " components", sep=""))
      mixtures[[k]] <- tryCatch({
        normalmixEM(tics[train,i], k=k, maxit=400,epsilon=1e-2,arbvar=FALSE)
      }, error = function(e) {
        return(-10000.0)
      }, finally = {
      })
        bics[k] <- BIC.mix(mixtures[[k]])
      #loglikes[k] <- loglike.normalmix(tics[test,i],mixture=mixture)
      #training.loglikes[k] <- loglike.normalmix(tics[train,i], mixture=mixture)
    }
    pdf(paste(dir, "/tIC_", i, "_BIC.pdf", sep=""))
    plot(bics)
    dev.off()
    #print(loglikes)
    #print(training.loglikes)
    print("finished computing training models")
    
    n.comp <- which.min(bics)
    #n.comp <- which.min(-1.0*loglikes)
    print(paste("choosing model with: ", n.comp, "components"))
    best.components[i] <- n.comp
    
    print("fitting model with that many number of components")
    best.model <- mixtures[[n.comp]]
    print("fit model, plotting now...")
    plot.mixture(tics[,i,drop=F], best.model, n.comp, filename = paste(dir, "/tIC", i, "_mixture_model.pdf", sep=""))
    
    if(n.comp > 1) {
        tic.classes <- apply(best.model$posterior, 1, which.max)
    } else {
        tic.classes <- as.factor(rep(1,dim(tics)[1]))
    }

    tica.classes[,i] <- tic.classes
  }
  }
  write.csv(tica.classes, paste(dir, "/tica_classes.csv", sep=""), sep=",")


}





compute.tic.splits <- function(tica, save.dir) {
  for(j in 1:dim(tica)[2]) {

  }
}

compute.decision.trees <- function(tica, tica.csv, tica.classes.csv, js, feature.npys, feature.csv = "", feature.residues.csv, save.dir) {
    print("js:")
    print(js)
    print("Starting")
     if(feature.csv != "") {
      print(feature.csv)
      features = my.read.table(feature.csv)
       
    } else {
      print(feature.csv)
      features = load.mat.matrix(feature.npys)
    }
    print("loaded features")
    feature.names <- get.feature.names(feature.residues.csv)
    colnames(features) <- feature.names
    print("added feature names to features")

    if(tica.csv != "") {
      tica = read.csv(tica.csv, stringsAsFactors = F, header = T, row.names = 1)
    } else {
      print("given tica")
    }
    print("read in tica coords")
    
    tica.classes <- read.csv(file = tica.classes.csv, sep = ",", header = T, row.names=1)
    
    tica.classes <- as.data.frame(unclass(tica.classes))
    tica.classes[] <- lapply(tica.classes,factor)
    print("read tica classes")

    foreach(k=1:length(js)) %dopar% {
      j = js[k]
      tica.j = tica[,j]
      tica.j.classes = tica.classes[,j]
      splits = find.class.divisions(tica.j, tica.j.classes)
      print("Found class divisions")
      new.classes <- rep(NA,length(tica.j))     
      splits <- sapply(splits, round, 1)
      for(i in 1:length(splits)) {
        if(i == 1) {
          rows = which(tica.j <= splits[i])
            print(rows[1:10])
          new.classes[rows] = paste("t<", splits[i],sep="")
        } else if(i==length(splits)) {
          rows = which(tica.j >= splits[i])
          print(rows[1:10])
          new.classes[rows] = paste(">t", splits[i],sep="")
          rows = which((tica.j >= splits[i-1]) & (tica.j <= splits[i]))
           print(rows[1:10])
          new.classes[rows] = paste(splits[i-1], ",", splits[i],sep="")
        } else {
          rows = which((tica.j >= splits[i-1]) & (tica.j <= splits[i]))
            print(rows[1:10])
          new.classes[rows] = paste(splits[i-1], ",", splits[i],sep="")
        }
      }
      new.classes <- factor(new.classes)
      tica.classes = new.classes
      #tica.classes = tica.j.classes
      features.tica <- as.data.frame(cbind(features,tica.classes))
      features.tica[,"tica.classes"] <- factor(features.tica[,"tica.classes"])
      print("combined data frames")
      #tica.classes = tica.j.classes
      #colnames(features.tica) <- make.names(colnames(features.tica))
      #print("combined data frames")
      #ticnums <- seq(1,dim(tica)[2],1)
      #ticnames <- sapply(ticnums, function(t) return(paste("tIC.", t, sep="")))
      #print(ticnames[1:10])
      #colnames(tica) <- ticnames
      
      #train.rows = sample(1:length(tica),length(tica)/2,replace=F)
      #train = features.tica[train.rows,]
      #test = features.tica[-train.rows,]
      #xtrain = features[train.rows,]
      #ytrain = tic.classes[train.rows]
      #xtest = features[-train.rows,]
      #ytest = tic.classes[-train.rows]

      tree.tica = tree(tica.classes ~., data=features.tica)
      pdf(paste(save.dir, "/decision_tree_initial_tIC", j, ".pdf", sep=""),width=10,height=5)
      plot(tree.tica, main = paste("Feature Decision Tree, tIC.", j, sep=""))
      text(tree.tica,pretty=1)
      dev.off()
      print("Fit first tree for")
      print(j)

      cv.tica = cv.tree(tree.tica, FUN = prune.misclass)
      print("completed cross-validation")
      pdf(paste(save.dir, "/cv_tree_tIC", j, ".pdf", sep=""))
      plot(cv.tica)
      dev.off()

      print("plotting best tree by cross-validation")
      num.nodes = cv.tica$size[which.min(cv.tica$dev)]
      prune.tica = prune.misclass(tree.tica,  best = num.nodes)
      pdf(paste(save.dir, "/decision_tree_pruned_tIC", j, ".pdf", sep=""), width = 10, height = 5)
      plot(prune.tica)
      text(tree.tica,pretty=1)
      dev.off()

      print("Fit tree")
  }
}

compute.rf <- function(tica, tica.csv, tica.classes.csv, js, feature.npys, feature.csv = "", feature.residues.csv, save.dir) {
    print("js:")
    print(js)
    print("Starting")
     if(feature.csv != "") {
      print(feature.csv)
      features = my.read.table(feature.csv)
       
    } else {
      print(feature.csv)
      features = load.mat.matrix(feature.npys)
    }
    print("loaded features")
    feature.names <- get.feature.names(feature.residues.csv)
    colnames(features) <- feature.names
    print("added feature names to features")

    if(tica.csv != "") {
      tica = read.csv(tica.csv, stringsAsFactors = F, header = T, row.names = 1)
    } else {
      print("given tica")
    }
    print("read in tica coords")
    
    tica.classes <- read.csv(file = tica.classes.csv, sep = ",", header = T, row.names=1)
    
    tica.classes <- as.data.frame(unclass(tica.classes))
    tica.classes[] <- lapply(tica.classes,factor)
    print("read tica classes")

    foreach(k=1:length(js)) %dopar% {
      j = js[k]
      tica.j = tica[,j]
      tica.j.classes = tica.classes[,j]
      splits = find.class.divisions(tica.j, tica.j.classes)
      print("Found class divisions")
      new.classes <- rep(NA,length(tica.j))     
      splits <- sapply(splits, round, 1)
      for(i in 1:length(splits)) {
        if(i == 1) {
          rows = which(tica.j <= splits[i])
            print(rows[1:10])
          new.classes[rows] = paste("t<", splits[i],sep="")
        } else if(i==length(splits)) {
          rows = which(tica.j >= splits[i])
          print(rows[1:10])
          new.classes[rows] = paste(">t", splits[i],sep="")
          rows = which((tica.j >= splits[i-1]) & (tica.j <= splits[i]))
           print(rows[1:10])
          new.classes[rows] = paste(splits[i-1], ",", splits[i],sep="")
        } else {
          rows = which((tica.j >= splits[i-1]) & (tica.j <= splits[i]))
            print(rows[1:10])
          new.classes[rows] = paste(splits[i-1], ",", splits[i],sep="")
        }
      }
      new.classes <- factor(new.classes)
      tica.classes = new.classes
      #tica.classes = tica.j.classes
      features.tica <- as.data.frame(cbind(features,tica.classes))
      features.tica[,"tica.classes"] <- factor(features.tica[,"tica.classes"])
      print("combined data frames")
      #tica.classes = tica.j.classes
      #colnames(features.tica) <- make.names(colnames(features.tica))
      #print("combined data frames")
      #ticnums <- seq(1,dim(tica)[2],1)
      #ticnames <- sapply(ticnums, function(t) return(paste("tIC.", t, sep="")))
      #print(ticnames[1:10])
      #colnames(tica) <- ticnames
      
      #train.rows = sample(1:length(tica),length(tica)/2,replace=F)
      #train = features.tica[train.rows,]
      #test = features.tica[-train.rows,]
      #xtrain = features[train.rows,]
      #ytrain = tic.classes[train.rows]
      #xtest = features[-train.rows,]
      #ytest = tic.classes[-train.rows]

      set.seed(2015)
      rf.tica = randomForest(features.tica[,-c(dim(features.tica)[2])],features.tica[,"tica.classes"],ntree=500,mtry=floor(sqrt(dim(features.tica)[2])),do.trace=10)
      pdf(paste(save.dir, "/tIC", j, "rf_importance.pdf", sep=""))
      varImpPlot(rf.tica)
      dev.off()

      save(rf.tica, file=paste(save.dir,"/tIC",j,"rf.rda",sep=""))

      pdf(paste(save.dir, "/tIC", j, "rf_oob_error.pdf", sep=""))
      plot(rf.tica$err.rate)
      dev.off()


      print("Fit and saved forest")
  }
}

write.rf.gini <- function(filename, new_filename) {
  load(filename)
  write.csv(rf.tica$importance, file = new_filename)
}

sample.tic.ranges <- function(tica.coords.csv, tica.classes.csv, js, samples.csv) {
    print("js:")
    print(js)
    print("Starting")

    if(tica.coords.csv != "") {
      tica = read.csv(tica.coords.csv, stringsAsFactors = F, header = T, row.names = 1)
    } else {
      print("given tica")
    }
    print("read in tica coords")
    
    tica.classes <- read.csv(file = tica.classes.csv, sep = ",", header = T, row.names=1)
    
    tica.classes <- as.data.frame(unclass(tica.classes))
    tica.classes[] <- lapply(tica.classes,factor)
    print("read tica classes")

    classes.per.tic = list()

    for(k in 1:length(js)) {
      j = js[k]
      tica.j = tica[,j]
      tica.j.classes = tica.classes[,j]
      splits = find.class.divisions(tica.j, tica.j.classes)
      print("Found class divisions")
      new.classes <- rep(NA,length(tica.j))     
      splits <- sapply(splits, round, 1)
      classes.per.tic.j  = rep(NA,length(splits)+1)
      for(i in 1:length(splits)) {
        if(i == 1) {
          rows = which(tica.j <= splits[i])
            print(rows[1:10])
          class.name = paste("t<", splits[i],sep="")   
          new.classes[rows] = class.name
          classes.per.tic.j[i] = class.name         
        } else if(i==length(splits)) {
          rows = which(tica.j >= splits[i])
          print(rows[1:10])
          class.name = paste("t>", splits[i],sep="")
          new.classes[rows] = class.name
          classes.per.tic.j[i+1] = class.name
          rows = which((tica.j >= splits[i-1]) & (tica.j <= splits[i]))
           print(rows[1:10])
           class.name = paste(splits[i-1], ",", splits[i],sep="")
          new.classes[rows] = class.name
          classes.per.tic.j[i] = class.name
        } else {
          rows = which((tica.j >= splits[i-1]) & (tica.j <= splits[i]))
            print(rows[1:10])
            class.name = paste(splits[i-1], ",", splits[i],sep="")
          new.classes[rows] = class.name
          classes.per.tic.j[i] = class.name
        }
        classes.per.tic[[j]] = classes.per.tic.j
        tica.classes[,j] <- new.classes
      }
    }

    print(classes.per.tic)
    tIC2.lowest.class <- classes.per.tic[[2]][length(classes.per.tic[[2]])]
    tIC2.medium.class <- classes.per.tic[[2]][length(classes.per.tic[[2]])-1]
    tIC2.higher.class <- classes.per.tic[[2]][length(classes.per.tic[[2]])-2]
    tIC7.highest.class <- classes.per.tic[[7]][length(classes.per.tic[[7]])]
    tIC7.higher.class <- classes.per.tic[[7]][length(classes.per.tic[[7]])-1]
    tIC9.highest.class <- classes.per.tic[[9]][length(classes.per.tic[[9]])]
    tIC9.higher.class <- classes.per.tic[[9]][length(classes.per.tic[[9]])-1]

    samples <- as.data.frame(cbind(rep(NA,25), rep(NA,25)))
    rows = which(tica.classes[,2] == tIC2.lowest.class & tica.classes[,7] == tIC7.higher.class & tica.classes[,9] == tIC9.higher.class)
    sample.rows <- sample(rows,5,replace=F)
    samples[seq(1,5,1),1] <- rep(paste("tIC2,",tIC2.lowest.class),5)
    samples[seq(1,5,1),2] <- rownames(tica)[sample.rows]
    rows = which(tica.classes[,2] == tIC2.medium.class & tica.classes[,7] == tIC7.higher.class & tica.classes[,9] == tIC9.higher.class)
    sample.rows <- sample(rows,5,replace=F)
    samples[seq(6,10,1),1] <- rep(paste("tIC2,",tIC2.medium.class),5)
    samples[seq(6,10,1),2] <- rownames(tica)[sample.rows]
    rows = which(tica.classes[,2] == tIC2.higher.class & tica.classes[,7] == tIC7.higher.class & tica.classes[,9] == tIC9.higher.class)
    sample.rows <- sample(rows,5,replace=F)
    samples[seq(11,15,1),1] <- rep(paste("tIC2,",tIC2.higher.class),5)
    samples[seq(11,15,1),2] <- rownames(tica)[sample.rows]
    rows = which(tica.classes[,2] == tIC2.higher.class & tica.classes[,7] == tIC7.highest.class & tica.classes[,9] == tIC9.higher.class)
    sample.rows <- sample(rows,5,replace=F)
    samples[seq(16,20,1),1] <- rep(paste("tIC7,",tIC7.highest.class),5)
    samples[seq(16,20,1),2] <- rownames(tica)[sample.rows]
    rows = which(tica.classes[,2] == tIC2.higher.class & tica.classes[,7] == tIC7.higher.class & tica.classes[,9] == tIC9.highest.class)
    sample.rows <- sample(rows,5,replace=F)
    samples[seq(21,25,1),1] <- rep(paste("tIC9,",tIC9.highest.class),5)
    samples[seq(21,25,1),2] <- rownames(tica)[sample.rows]
    write.csv(samples,samples.csv)
    rows = which(tica.classes[,2] == tIC2.higher.class & tica.classes[,7] == tIC7.highest.class & tica.classes[,9] == tIC9.highest.class)
    sample.rows <- sample(rows,5,replace=F)
    samples[seq(26,30,1),1] <- rep(paste("tIC9,",tIC9.highest.class,"tIC7,",tIC7.highest.class),5)
    samples[seq(26,30,1),2] <- rownames(tica)[sample.rows]
    write.csv(samples,samples.csv)

  
}


analyze.extreme.tic.values <- function(low.csv, high.csv, feature.residues.csv, i, dir) {
    low.values <- read.csv(low.csv, stringsAsFactors = F, header = F, colClasses = "numeric")
    #print(low.values[1:10,1:10])
    high.values <- read.csv(high.csv, stringsAsFactors = F, header = F, colClasses = "numeric")
   #print(feature.residues.csv)
    feature.residues <- 
    #print(feature.residues[1:10,])
    feature.names <- combine.columns(feature.residues[,2],feature.residues[,3])
    colnames(low.values) <- feature.names
    colnames(high.values) <- feature.names 
    low.means <- apply(low.values,2,mean)
    high.means <- apply(high.values,2,mean)
    differences <- abs(high.means-low.means)
   #print(differences[1:50])
    differences <- differences[order(-1.0*differences)]
   #print(differences[1:50])
    #barplot(differences[1:50],las=2, main = paste("tIC",i, " feature coefficient", sep="")) 
    pdf(paste(dir, "/greatest_feature_differences_tIC.", i, ".pdf", sep=""),width=20, height=5)
    barplot(differences[1:50],las=2, main = paste("tIC",i, " top 50 feature changes", sep=""))
    dev.off()
    
    if(1 == 1) {
    low.values.1 <- low.values
    low.values.2 <- low.values
    colnames(low.values.1) <- feature.residues[,2]
    colnames(low.values.2) <- feature.residues[,3]
    low.values <- rbind(t(low.values.1), t(low.values.2))
    
    high.values.1 <- high.values
    high.values.2 <- high.values
    colnames(high.values.1) <- feature.residues[,2]
    colnames(high.values.2) <- feature.residues[,3]
    high.values <- rbind(t(high.values.1), t(high.values.2))
    
    low.averages <- cluster.averages.new(low.values)
    high.averages <- cluster.averages.new(high.values)
    
    low.means <- apply(low.averages, 1, mean)
    high.means <- apply(high.averages, 1, mean)
    differences <- abs(high.means-low.means)
    differences <- differences[order(-1.0*differences)]
    
    pdf(paste(dir, "/greatest_residue_differences_tIC.", i, ".pdf", sep=""), width=20, height=5)
    barplot(differences[1:50],las=2, main = paste("tIC",i, " greatest residue difference", sep=""))
    dev.off()
    }
}

analyze.tic.feature.correlations <- function(info.csv = "", feature.residues.csv = "", dir, file_string, title_string) {
  info = read.csv(info.csv, stringsAsFactors = F, header = F, colClasses = "numeric")
  names <- sapply(seq(1,dim(info)[2]), function(x) paste("tIC.", x, sep=""))
  colnames(info) <- names
  feature.residues <- data.frame(read.csv(feature.residues.csv, stringsAsFactors = F, header=T))
  feature.names <- combine.columns(feature.residues[,2],feature.residues[,3])
 #print(feature.names[1:10])
  rownames(info) <- feature.names
  for(i in 1:dim(info)[2]) {
   #print(i)
    tic <- names[i]
   #print(tic)
    infos <- info[,i,drop=F]
    infos <- infos[order(-1.0*infos[,1]),1,drop=F]
    infos.barplot <- infos[,1]
    names(infos.barplot) <- rownames(infos)
    
    pdf(paste(dir, "/", tic, file_string, ".pdf", sep=""), width=20, height=5)
    barplot(infos.barplot[1:50], las=2, main = paste(tic, " ", title_string, sep=""))
    dev.off()
  }
  
  if(1 == 2) {
  info.1 <- info
  info.2 <- info
  rownames(info.1) <- feature.residues[,2]
  rownames(info.2) <- feature.residues[,3]
  helix <- which.helix(feature.residues[,2], only.helices=TRUE)
 #print("Helix:")
 #print(dim(helix))
 #print(helix[1:10,])
  helix.averages <- cluster.averages.new(helix)
  helix.sorted <- sort.by.column(helix.averages)
  for(i in 1:dim(tic.residue.averages)[2]) {
    pdf(paste(analysis.dir, "/barplot_average_helix_coefficient_tIC", i, ".pdf", sep=""))
    barplot(helix.sorted[[i]][,i],las=2, main = paste("tIC",i, " average helix coefficient", sep=""))
    dev.off()
  }
  
  tic.residues <- get.tic.residues.pairs(analysis.dir = analysis.dir, features.dir = features.dir, tica = tica)
  interhelix <- interhelix.distance(tic.residues)
  interhelix.averages <- cluster.averages.new(interhelix)
  interhelix.averages.sorted <- sort.by.column(interhelix.averages)
  for(i in 1:dim(interhelix.averages)[2]) {
    pdf(paste(analysis.dir, "/barplot_average_inter_helix_coefficient", i, ".pdf", sep=""))
    barplot(interhelix.averages.sorted[[i]][,i],las=2, main = paste("tIC",i, " average interhelix distance coefficient", sep=""))
    dev.off()
  }
  }
  
}



#docking.vs.tica.lm.cv <- cv.glm(docking.tica, docking.vs.tica.lm,K=10)
#summary(docking.vs.tica.lm)
#docking.vs.tica.lm.cv
#docking.vs.tica.lm.cv$delta
#docking.binary.tica <- docking.tica
#replacement <- rep(0,dim(docking.binary.tica)[1])
#replacement[docking.binary.tica[,1] > (mean(docking.binary.tica[,1])+1.0*sd(docking.binary.tica[,1]))] <- 1.0
#docking.binary.tica[,1] <- replacement
#docking.binary.vs.tica.glm <- glm(rel.docking.score~.,data=docking.binary.tica,family=binomial)
#docking.binary.vs.tica.glm.cv <- cv.glm(docking.binary.tica, docking.binary.vs.tica.glm,K=5)
#summary(docking.vs.tica.glm)
#docking.binary.vs.tica.glm.cv
#docking.binary.tica <- docking.binary.tica[order(-1.0*docking.binary.tica[,10]),]
#docking.binary.vs.tica9.glm <- glm(rel.docking.score~tIC.9,data=docking.binary.tica,family=binomial)
#docking.binary.vs.tica9.glm
#predicted.glm <- predict(docking.binary.vs.tica9.glm, docking.binary.tica[,10,drop=F], type="response")
#plot(docking.binary.tica[,"tIC.9"],docking.binary.tica[,1])
#lines(docking.binary.tica[,"tIC.9"], predicted.glm)
#docking.binary.vs.tica.lm <- lm(rel.docking.score~.,data=docking.binary.tica)
#summary(docking.binary.vs.tica.lm)
#pairs(docking.tica)
#