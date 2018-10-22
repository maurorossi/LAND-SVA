#########################################################################
#########################################################################
####                                                                 ####
####                                                                 ####
####                                                                 ####
####                           LAND-SVA                              ####
####          LANDSLIDE SUSCEPTIBILITY VARIABLE ANALYSIS             ####
####                           IRPI CNR                              ####
####                    MAURO ROSSI - IRPI CNR                       ####
####                  TXOMIN BORNAETXEA - UPV/EHU                    ####
####                                                                 ####
####                     v1r0b1 - 11 November 2016                   ####
#### Copyright (C) 2016 Mauro Rossi, Txomin Bornaetxea               ####
####                                                                 ####
#### This program is free software; you can redistribute it and/or   ####
#### modify it under the terms of the GNU General Public License     ####
#### as published by the Free Software Foundation; either version 2  ####
#### of the License, or (at your option) any later version. ###      ####
####                                                                 ####
#### This program is distributed in the hope that it will be useful, ####
#### but WITHOUT ANY WARRANTY; without even the implied warranty of  ####
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the ####
#### GNU General Public License for more details.                    ####
####                                                                 ####
####     Istituto di Ricerca per la Protezione Idrogeologica         ####
####              Consiglio Nazionale delle Ricerche                 ####
####                    Gruppo di Geomorfologia                      ####
####                  Via della Madonna Alta, 126                    ####
####                    06128 Perugia (Italia)                       ####
####                       +39 075 5014421                           ####
####                       +39 075 5014420                           ####
####                   mauro.rossi@irpi.cnr.it                       ####
####                  geomorfologia@irpi.cnr.it                      ####
####                                                                 ####
####    Universidad del País Vasco/Euskal Herriko Unibertsitatea     ####
####       Facultad de Ciencias, Departamento de Geodinamica         ####
####               Zientzia Fakultatea, Geodinamika Saila            ####
####                     txomin.bornaetxea@ehu.eus                   ####
####                                                                 ####
####           This script was prepared using R 3.1.3                ####
####         The script requires the following R packages:           ####
####                       1: corrplot                               ####
####                       2: perturb                                ####
####                       3: Hmisc                                  ####
####                       4: data.table                             ####
####                       5: RColorBrewer                           ####
####                       6: rgdal                                  ####
####                       7:                                        ####
####                       8:                                        ####
####                       9:                                        ####
####                                                                 ####
####     INPUTS: 1) datatable_inventory.RData                        ####
####                produced by the script LAND-SIP.R                ####

####                                                                 ####
####                                                                 ####
#########################################################################
#########################################################################

#R CMD BATCH --no-save --no-restore '--args -wd /media/disco_dati/R/grid_to_xyz/generali/adb_07_medium/' LAND-SVA_v1r0b1_20161111.R variable_analysis.log

rm(list=(ls()))
graphics.off()
#setwd("X:/R/grid_to_xyz/")
pars <-commandArgs(trailingOnly=TRUE)

if (length(table(pars == "-wd"))==2)
  {
  wd_selected<-pars[which(pars=="-wd")+1]
  } else
  {
  wd_selected<-"X:/R/grid_to_xyz/LAND-SVA_github/" # manually specified 
  }
setwd(wd_selected)
#memory.limit(size=12000)

### --------------- SVA analysis paramter definition --------------- ###

rdata_file<-paste("datatable_inventory.RData",sep="")
load(file=rdata_file)

enable_NA_removal<-TRUE # If TRUE rows with at least an NA value will be removed. If FALSE and error message will be returned.
enable_multicollinearity_test<-TRUE
type_correlation<-"pearson" # It could be "pearson" or "spearman" and it specifies the type of correlations to compute. Spearman correlations are the Pearson linear correlations computed on the ranks of non-missing elements, using midranks for ties. Pearson's coefficient and Spearman's rank order coefficient each measure aspects of the relationship between two variables. They are closely related, but not the same. Spearman's coefficient measures the rank order of the points. Pearson's coefficient measures the linear relationship between the two.

export_shapefiles<-TRUE # Enable this to export shapefile of points corresponding to the data tables, Usefull to check location of point of training and validation datasets using GIS clients
export_txtfiles<-TRUE # Enable this to export the training and validation tables in tab separeted .txt format


### --------------- Data conversion --------------- ###

library(data.table)

# converting to data.table format
training.table<-data.table(training.table)
original_rows_training<-dim(training.table)[1]
validation.table<-data.table(validation.table)
original_rows_validation<-dim(validation.table)[1]


### --------------- NAs analysis and selection --------------- ###

index_selection_training<-is.finite(rowSums(training.table))
training.table<-training.table[index_selection_training,1:dim(training.table)[2],with=FALSE]
variables_training<-training.table[,3:dim(training.table)[2],with=FALSE]

index_selection_validation<-is.finite(rowSums(validation.table))
validation.table<-validation.table[index_selection_validation,1:dim(validation.table)[2],with=FALSE]
variables_validation<-validation.table[,3:dim(validation.table)[2],with=FALSE]


### ---------------  Conditional Density Plots --------------- ###
for(count in 1:dim(variables_training)[2])
  {
  #count<-2
  selected_variable<-names(variables_training)[count]
  print(paste("Conditional and Spine plots training variable: ",selected_variable," - Done: ",round(count/dim(variables_training)[2]*100,1),"%",sep=""))
  pdf(paste("ConditionalDensityPlot_",selected_variable,".pdf",sep=""))
  cdplot(x=as.numeric(variables_training[,count,with=FALSE][[1]]),y=as.factor(training.table[,2,with=FALSE][[1]]),bw="nrd0",kernel="gaussian",xlab=selected_variable,ylab="Dependent variable",main=paste("Conditional plot: ",selected_variable,sep=""))
  dev.off()
  
  pdf(paste("SpinePlot_",selected_variable,".pdf",sep=""))
  breaks_sel<-nclass.Sturges(as.numeric(variables_training[,count,with=FALSE][[1]]))
  #breaks_sel<-nclass.scott(as.numeric(variables_training[,count,with=FALSE][[1]]))
  #breaks_sel<-nclass.FD(as.numeric(variables_training[,count,with=FALSE][[1]]))
  spineplot(x=as.numeric(variables_training[,count,with=FALSE][[1]]),y=as.factor(training.table[,2,with=FALSE][[1]]),breaks=breaks_sel,xlab=selected_variable,ylab="Dependent variable",main=paste("Spineplot: ",selected_variable,sep=""))
  dev.off()
  
  pdf(paste("DensityPlot_",selected_variable,".pdf",sep=""))
  dependent_values<-as.numeric(names(table(training.table[,2,with=FALSE])))
  density_results_names<-paste("den_",dependent_values,sep="")
  require(RColorBrewer)
  #display.brewer.all()
  colors_vector<-brewer.pal(length(dependent_values)+1,"Set1")
  range_y_den<-c(0,0.01)
  for(count_den in 1:length(dependent_values))
    {
    #count_den<-1
    dependent_values_selected<-dependent_values[count_den]
    index_dependent<-which(training.table[,2,with=FALSE]==dependent_values_selected)
    assign(density_results_names[count_den],density(x=as.numeric(variables_training[index_dependent,count,with=FALSE][[1]]),bw="nrd0",kernel="gaussian"))
    if(max(get(density_results_names[count_den])$y,na.rm=TRUE)>range_y_den[2]) range_y_den<-c(0,max(get(density_results_names[count_den])$y,na.rm=TRUE))
    }
  plot(NULL,NULL,xlab=selected_variable,ylab="Density",xlim=range(as.numeric(variables_training[,count,with=FALSE][[1]])),ylim=range_y_den,main=paste("Density plot: ",selected_variable,sep=""))
  for(count_plot_den in 1:length(dependent_values))
    {
    lines(get(density_results_names[count_plot_den]), col=colors_vector[count_plot_den],lwd=2)
    }
  legend("topleft",legend=dependent_values,lty=1,lwd=2,col=colors_vector,bty="n") #
  dev.off()
  }



### --------------- Multicolinearity test for the training and validation dataset --------------- ###

if(enable_multicollinearity_test==TRUE) 
  {
  #load collinearity package (perturb)
  library(perturb)
  #colnames(training.table)
  collinearity.test.training<-colldiag(variables_training)
  collinearity.test.validation<-colldiag(variables_validation)
  #collinearity.test.training$condindx 
  #collinearity.test.training$pi 
  #range(collinearity.test.training$condindx)
  
  if(range(collinearity.test.training$condindx)[2] >= 30)
    {      
    collinearity.value.training<-"Some explanatory variables are collinear"
    } else {      
    collinearity.value.training<-"Explanatory variables are not collinear"
    }
  
  print(collinearity.test.training,fuzz=.5)
  collinearity.evaluation.matrix.training<-print(collinearity.test.training,fuzz=.5)
  write.table("COLLINEARITY ANALYSIS RESULT",file="Variables_CollinearityAnalysis_training.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("EXPLANATION",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("This analysis was performed with Colldiag an implementation of the regression collinearity diagnostic procedures found in Belsley, Kuh,
              and Welsch (1980). These procedures examine the ?conditioning? of the matrix of independent variables. The procedure computes the condition
              indexes of the matrix. If the largest condition index (the condition number) is large (Belsley et al suggest 30 or higher), then there may be
              collinearity problems. All large condition indexes may be worth investigating. The procedure also provides further information that may help to
              identify the source of these problems, the variance decomposition proportions associated with each condition index. If a large condition
              index (> 30) is associated with two or more variables with large variance decomposition proportions, these variables may be causing collinearity problems.
              Belsley et al suggest that a large proportion is 50 percent or more.",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t",
              row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("RESULTS",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(paste("Largest condition index (the condition number) =",range(collinearity.test.training$condindx)[2]),file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(collinearity.value.training,file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Matrix of the variance decomposition proportions associated with each condition index (1st column)",file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(colnames(collinearity.evaluation.matrix.training),collinearity.evaluation.matrix.training),file="Variables_CollinearityAnalysis_training.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  
  
  if(range(collinearity.test.validation$condindx)[2] >= 30)
    {      
    collinearity.value.validation<-"Some explanatory variables are collinear"
    } else {      
    collinearity.value.validation<-"Explanatory variables are not collinear"
    }
  
  print(collinearity.test.validation,fuzz=.5)
  collinearity.evaluation.matrix.validation<-print(collinearity.test.validation,fuzz=.5)
  write.table("COLLINEARITY ANALYSIS RESULT",file="Variables_CollinearityAnalysis_validation.txt", quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("EXPLANATION",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("This analysis was performed with Colldiag an implementation of the regression collinearity diagnostic procedures found in Belsley, Kuh,
              and Welsch (1980). These procedures examine the ?conditioning? of the matrix of independent variables. The procedure computes the condition
              indexes of the matrix. If the largest condition index (the condition number) is large (Belsley et al suggest 30 or higher), then there may be
              collinearity problems. All large condition indexes may be worth investigating. The procedure also provides further information that may help to
              identify the source of these problems, the variance decomposition proportions associated with each condition index. If a large condition
              index (> 30) is associated with two or more variables with large variance decomposition proportions, these variables may be causing collinearity problems.
              Belsley et al suggest that a large proportion is 50 percent or more.",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t",
              row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("RESULTS",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(paste("Largest condition index (the condition number) =",range(collinearity.test.validation$condindx)[2]),file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(collinearity.value.validation,file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table("Matrix of the variance decomposition proportions associated with each condition index (1st column)",file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  write.table(rbind(colnames(collinearity.evaluation.matrix.validation),collinearity.evaluation.matrix.validation),file="Variables_CollinearityAnalysis_validation.txt", append=TRUE, quote = FALSE,sep = "\t", row.names=FALSE, col.names=FALSE)
  
  }


### --------------- Calculating and plotting correlation matrix --------------- ###
library(Hmisc)
coeff_training<-rcorr(as.matrix(variables_training),type=type_correlation) # This comand run the correlation matrix and the P value matrix with the Pearson method and exculding the NA values
coeff_validation<-rcorr(as.matrix(variables_validation),type=type_correlation) # This comand run the correlation matrix and the P value matrix with the Pearson method and exculding the NA values


library(corrplot)
colors<-colorRampPalette(c("darkred","grey40","forestgreen"))(100)

pdf("Variables_Correlogram_matrix_training.pdf") 
corrplot.mixed(coeff_training$r, lower="ellipse", upper="number",tl.col="black",tl.pos="lt",number.cex=0.7,tl.cex=0.6,cl.cex=0.6,cl.ratio=0.2,cl.align.text = "l",title = "Correlation matrix", order = "original",p.mat=coeff_training$P,sig.level=0.01,insig="blank")
#corrplot.mixed(coeff_training$r, lower="ellipse", upper="number",tl.col="black",tl.pos="lt",col=colors,number.cex=0.7,tl.cex=0.6,cl.cex=0.6,cl.ratio=0.2,cl.align.text = "l",title = "Correlation matrix", order = "original",p.mat=coeff_training$P,sig.level=0.01,insig="blank")
dev.off()

pdf("Variables_Correlogram_matrix_validation.pdf") 
corrplot.mixed(coeff_validation$r, lower="ellipse", upper="number",tl.col="black",tl.pos="lt",number.cex=0.7,tl.cex=0.6,cl.cex=0.6,cl.ratio=0.2,cl.align.text = "l",title = "Correlation matrix", order = "original",p.mat=coeff_validation$P,sig.level=0.01,insig="blank")
#corrplot.mixed(coeff_validation$r, lower="ellipse", upper="number",tl.col="black",tl.pos="lt",col=colors,number.cex=0.7,tl.cex=0.6,cl.cex=0.6,cl.ratio=0.2,cl.align.text = "l",title = "Correlation matrix", order = "original",p.mat=coeff_validation$P,sig.level=0.01,insig="blank")
dev.off()


### --------------- NAs removal and writing output training and validation tables --------------- ###

if(enable_NA_removal==TRUE)
  {
  if(original_rows_training!=dim(training.table)[1]) print("Warning: Analysis will be executed excluding rows with NA from the training set")
  if(original_rows_validation!=dim(validation.table)[1]) print("Warning: Analysis will be executed excluding rows with NA from the validation set")
  
  training.table<-data.frame(training.table)
  validation.table<-data.frame(validation.table)

  training.xy.table@data<-data.frame(id=training.xy.table@data[index_selection_training,])
  training.xy.table@coords<-training.xy.table@coords[index_selection_training,]
  
  validation.xy.table@data<-data.frame(id=validation.xy.table@data[index_selection_validation,])
  validation.xy.table@coords<-validation.xy.table@coords[index_selection_validation,]
  
  save(list=c("training.table","validation.table","training.xy.table","validation.xy.table"),file = paste("datatable_inventory.RData",sep=""))
  
  } else
  {
  if(original_rows_training!=dim(training.table)[1]) print("Error: training set contains rows with NA values. All raster layers should contain finite values within the mask. If you want to execute the analysis excluding NA set the variable enable_NA_removal<-TRUE")
  if(original_rows_validation!=dim(validation.table)[1]) print("Error: validation set contains rows with NA values. All raster layers should contain finite values within the mask. If you want to execute the analysis excluding NA set the variable enable_NA_removal<-TRUE")
  }


result_dir_susceptibility<-getwd()
if(export_shapefiles==TRUE)
  {
  require(rgdal)
  writeOGR(training.xy.table,dsn=paste(result_dir_susceptibility,sep=""),layer="training",driver="ESRI Shapefile",overwrite_layer=TRUE)
  writeOGR(validation.xy.table,dsn=paste(result_dir_susceptibility,sep=""),layer="validation",driver="ESRI Shapefile",overwrite_layer=TRUE)
  }

if(export_txtfiles==TRUE)
  {
  write.table(training.table,file=paste(result_dir_susceptibility,"/training.txt",sep=""),sep="\t",dec=".",row.names=FALSE,col.names=TRUE)
  write.table(validation.table,file=paste(result_dir_susceptibility,"/validation.txt",sep=""),sep="\t",dec=".",row.names=FALSE,col.names=TRUE)
  }


