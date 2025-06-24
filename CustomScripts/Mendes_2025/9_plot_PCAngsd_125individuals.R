#plot the resulting PCA and Admixture plots after running PCAngsd and NGS Admix

########################## PCA map colours All individuals (all species) ############################
#125 individuals in angsd
#13.11.2023 new colours and symbols according to new map

#read the popmap
popmap<-read.table(file="popmap_125inds_lcWGS_scephalus_merged_PCA.txt", header = FALSE)
#get the names of the sampling locations 
pops<-unique(popmap[,2])
#get the number of sampling locations
npops<-as.numeric(length(unique(popmap[,2])))
#get the sampling location to which every individual belongs
pop_inds<-popmap[,2]

#read the cov matrix from PCAngsd
m<-as.matrix(read.table(file="merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.cov", header=FALSE, sep="", dec="."))
e<-eigen(m)
#create a "colour pallete"
c_carol_north<-"indianred2"
c_carol_mondego<-"lightpink"
c_lis<-"deepskyblue"
c_alcoa<-"lightpink"
c_pyr_tejo<-"deepskyblue"
c_zezere<-"cyan1"
c_nabao<-"deepskyblue"
c_lizandro<-"deepskyblue"
c_guadiana<-"blue2"
c_guadalquivir<-"navy"
c_guadalete<-"navy"
c_guadalhorce<-"grey50"
c_velez<-"grey50"
c_guadalfeo<-"navy"
c_sado<-"darkviolet"
c_smartinho<-"plum3"
c_torgalensis<-"yellow"
c_aradensis<-"green2"
c_malacitanus<-"chocolate4"
c_mijares<-"magenta"
c_jucar_pyr<-"deepskyblue"
c_jucar_val<-"magenta"
#define the colour of each individual
indcolours<-c(rep(c_carol_north,24),rep(c_carol_mondego,6),rep(c_lis,2),rep(c_alcoa,2),rep(c_pyr_tejo,14),rep(c_zezere,25),rep(c_nabao,2),rep(c_lizandro,2),rep(c_guadiana,10),rep(c_guadalquivir,4),rep(c_guadalete,2),rep(c_guadalhorce,2),rep(c_velez,2),rep(c_guadalfeo,2),rep(c_sado,4),rep(c_smartinho,9),rep(c_torgalensis,2),rep(c_aradensis,4),rep(c_malacitanus,2),rep(c_mijares,2),rep(c_jucar_pyr,1),rep(c_jucar_val,2))
#check we have colours for 125 individuals
length(indcolours)
#define the symbol of each individual
s_carol<-rep(0,30)
s_lis<-rep(1,2)
s_alcoa<-rep(0,2)
s_pyrnorte<-rep(1,43)
s_pyrsouth<-rep(2,22)
s_sado<-rep(5,13)
s_torgalensis<-rep(0,2)
s_aradensis<-rep(1,4)
s_malacitanus<-rep(5,2)
s_valentinus<-rep(2,5)
indsymbols<-c(s_carol,s_lis,s_alcoa,s_pyrnorte,s_pyrsouth,s_sado,s_torgalensis,s_aradensis,s_malacitanus,s_valentinus)
#check we have symbols for all individuals
length(indsymbols)
#define the colour of each population
popcolours<-c(c_carol_north,c_carol_mondego,c_lis,c_alcoa,c_pyr_tejo,c_zezere,c_nabao,c_lizandro,c_guadiana,c_guadalquivir,c_guadalete,c_guadalhorce,c_velez,c_guadalfeo,c_sado,c_smartinho,c_torgalensis,c_aradensis,c_malacitanus,c_mijares,c_jucar_pyr,c_jucar_val)
#check that we defined a colour for all populations, that is, popcolours has to be
#the same length of npops
length(popcolours)==npops
#define the symbols of each population
l_carol_north<-0
l_carol_mondego<-0
l_lis<-1
l_alcoa<-0
l_pyr_tejo<-1
l_zezere<-1
l_nabao<-1
l_lizandro<-1
l_guadiana<-2
l_guadalquivir<-2
l_guadalete<-2
l_guadalhorce<-2
l_velez<-2
l_guadalfeo<-2
l_sado<-5
l_smartinho<-5
l_torgalensis<-0
l_aradensis<-1
l_malacitanus<-5
l_mijares<-2
l_jucar_pyr<-2
l_jucar_val<-2
legsymbols<-c(l_carol_north,l_carol_mondego,l_lis,l_alcoa,l_pyr_tejo,l_zezere,l_nabao,l_lizandro,l_guadiana,l_guadalquivir,l_guadalete,l_guadalhorce,l_velez,l_guadalfeo,l_sado,l_smartinho,l_torgalensis,l_aradensis,l_malacitanus,l_mijares,l_jucar_pyr,l_jucar_val)
#confirm that we defined a symbol for each population
length(legsymbols)==npops
#print just the legend
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("bottomleft", pops, col=popcolours, pch=legsymbols, cex=1, ncol = 3)

#Create pdf to save the plots
pdf(file = "PCAngsd_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_ADMIXMAPCOLOURS_new_map_symbols_colours_13_11_2023.pdf", width = 11.693, height=8.268)
#print the PC1 and PC2 plot
plot(e$vectors[,1:2],ylab=paste0("PC2 (",round(e$values[2]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC1 (",round(e$values[1]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC1 and PC3 plot
plot(e$vectors[,c(1,3)],ylab=paste0("PC3 (",round(e$values[3]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC1 (",round(e$values[1]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC1 and PC4 plot
plot(e$vectors[,c(1,4)],ylab=paste0("PC4 (",round(e$values[4]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC1 (",round(e$values[1]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC2 and PC3 plot
plot(e$vectors[,2:3],ylab=paste0("PC3 (",round(e$values[3]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC2 (",round(e$values[2]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC2 and PC4 plot
plot(e$vectors[,c(2,4)],ylab=paste0("PC4 (",round(e$values[4]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC2 (",round(e$values[2]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC3 and PC4 plot
plot(e$vectors[,c(3,4)],ylab=paste0("PC4 (",round(e$values[4]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC3 (",round(e$values[3]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC1 and PC5 plot
plot(e$vectors[,c(1,5)],ylab=paste0("PC5 (",round(e$values[5]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC1 (",round(e$values[1]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC2 and PC5 plot
plot(e$vectors[,c(2,5)],ylab=paste0("PC5 (",round(e$values[5]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC2 (",round(e$values[2]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC3 and PC5 plot
plot(e$vectors[,c(3,5)],ylab=paste0("PC5 (",round(e$values[5]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC3 (",round(e$values[3]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)
#print the PC4 and PC5 plot
plot(e$vectors[,c(4,5)],ylab=paste0("PC5 (",round(e$values[5]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC4 (",round(e$values[4]/sum(e$values)*100,digits=2),"%)"),col=indcolours,pch=indsymbols,lwd=3,cex=2.5)

dev.off()

rm(list = ls())


##################### PCAngsd ADMIX 125inds All individuals (all species) ##############################

#PCAngsd Admix option was run for K values 2 to 17
###Plot the output of the various K values for 125 individuals

#create a pdf file to save the plots
pdf(file = "PCAngsd_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_K2_K17_ADMIX_REORDERED.pdf", width = 11.693, height=8.268)

#reorder the columns so the plots have the individuals ordered in a way that is easier to visualize the hybrids
angsd_order<-c(1:125)
paper_order<-c(c(1:44),c(49:75),c(45:48),c(76,77),c(104:112),c(100:103),c(78:99),c(113:125))


#K=2
temp_dataset_K2<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.2.Q"))
datasetK2<-t(temp_dataset_K2)
colnames(datasetK2)<-angsd_order
datasetK2_reordered<-datasetK2[,paper_order]
#barplot(dataset2,  col=c("grey", "black"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K2")
barplot(datasetK2_reordered,  col=c("grey", "black"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K2")

#K=3
temp_dataset_K3<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.3.Q"))
datasetK3<-t(temp_dataset_K3)
colnames(datasetK3)<-angsd_order
datasetK3_reordered<-datasetK3[,paper_order]
#barplot(dataset3,  col=c("black", "grey", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K3")
barplot(datasetK3_reordered,  col=c("black", "grey", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K3")

#K=4
temp_dataset_K4<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.4.Q"))
datasetK4<-t(temp_dataset_K4)
colnames(datasetK4)<-angsd_order
datasetK4_reordered<-datasetK4[,paper_order]
#barplot(datasetK4,  col=c("grey", "orangered", "deepskyblue", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K4")
barplot(datasetK4_reordered,  col=c("grey", "orangered", "deepskyblue", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K4")

#K=5
temp_dataset_K5<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.5.Q"))
datasetK5<-t(temp_dataset_K5)
colnames(datasetK5)<-angsd_order
datasetK5_reordered<-datasetK5[,paper_order]
#barplot(datasetK5,  col=c("grey", "lightgreen", "deepskyblue", "darkviolet", "orangered"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K5")
barplot(datasetK5_reordered,  col=c("grey", "lightgreen", "deepskyblue", "darkviolet", "orangered"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K5")

#K=6
temp_dataset_K6<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.6.Q"))
datasetK6<-t(temp_dataset_K6)
colnames(datasetK6)<-angsd_order
datasetK6_reordered<-datasetK6[,paper_order]
#barplot(datasetK6,  col=c("deepskyblue", "orangered", "darkviolet", "grey","yellow", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K6")
barplot(datasetK6_reordered,  col=c("deepskyblue", "orangered", "darkviolet", "grey","yellow", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K6")

#K=7
temp_dataset_K7<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.7.Q"))
datasetK7<-t(temp_dataset_K7)
colnames(datasetK7)<-angsd_order
datasetK7_reordered<-datasetK7[,paper_order]
#barplot(datasetK7,  col=c("darkviolet", "deepskyblue", "lightgreen", "lightpink","indianred2", "grey", "black"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K7")
barplot(datasetK7_reordered,  col=c("darkviolet", "deepskyblue", "lightgreen", "lightpink","indianred2", "grey", "black"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K7")

#K=8
temp_dataset_K8<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.8.Q"))
datasetK8<-t(temp_dataset_K8)
colnames(datasetK8)<-angsd_order
datasetK8_reordered<-datasetK8[,paper_order]
#barplot(datasetK8,  col=c("lightpink", "deepskyblue", "blue2", "indianred2", "navy", "lightgreen", "magenta", "darkviolet"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K8")
barplot(datasetK8_reordered,  col=c("lightpink", "deepskyblue", "blue2", "indianred2", "navy", "lightgreen", "magenta", "darkviolet"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K8")

#K=9
temp_dataset_K9<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.9.Q"))
datasetK9<-t(temp_dataset_K9)
colnames(datasetK9)<-angsd_order
datasetK9_reordered<-datasetK9[,paper_order]
#barplot(datasetK9,  col=c("blue2", "chocolate4", "darkviolet", "yellow", "indianred2", "lightpink", "magenta", "green2", "deepskyblue"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K9")
barplot(datasetK9_reordered,  col=c("blue2", "chocolate4", "darkviolet", "yellow", "indianred2", "lightpink", "magenta", "green2", "deepskyblue"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K9")

#K=10
temp_dataset_K10<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.10.Q"))
datasetK10<-t(temp_dataset_K10)
colnames(datasetK10)<-angsd_order
datasetK10_reordered<-datasetK10[,paper_order]
#barplot(datasetK10,  col=c("deepskyblue", "green2", "blue2", "indianred2","navy","chocolate4", "magenta", "lightpink", "darkviolet", "yellow"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K10")
barplot(datasetK10_reordered,  col=c("deepskyblue", "green2", "blue2", "indianred2","navy","chocolate4", "magenta", "lightpink", "darkviolet", "yellow"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K10")

#K=11
temp_dataset_K11<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.11.Q"))
datasetK11<-t(temp_dataset_K11)
colnames(datasetK11)<-angsd_order
datasetK11_reordered<-datasetK11[,paper_order]
#barplot(datasetK11,  col=c("deepskyblue", "green2", "navy", "magenta","yellow","lightpink", "darkviolet", "indianred2", "chocolate4", "blue2", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K11")
barplot(datasetK11_reordered,  col=c("deepskyblue", "green2", "navy", "magenta","yellow","lightpink", "darkviolet", "indianred2", "chocolate4", "blue2", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K11")

#K=12
temp_dataset_K12<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.12.Q"))
datasetK12<-t(temp_dataset_K12)
colnames(datasetK12)<-angsd_order
datasetK12_reordered<-datasetK12[,paper_order]
#barplot(datasetK12,  col=c("darkviolet", "lightpink", "indianred2", "navy","chocolate4", "yellow", "green2", "blue2", "deepskyblue", "cyan1", "grey50","magenta"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K12")
barplot(datasetK12_reordered,  col=c("darkviolet", "lightpink", "indianred2", "navy","chocolate4", "yellow", "green2", "blue2", "deepskyblue", "cyan1", "grey50","magenta"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K12")

#K=13
temp_dataset_K13<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.13.Q"))
datasetK13<-t(temp_dataset_K13)
colnames(datasetK13)<-angsd_order
datasetK13_reordered<-datasetK13[,paper_order]
#barplot(datasetK13,  col=c("cyan1", "lightpink", "deepskyblue","indianred2", "magenta", "green2", "blue2","grey50", "darkviolet","lightblue","chocolate4","navy", "yellow"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K13")
barplot(datasetK13_reordered,  col=c("cyan1", "lightpink", "deepskyblue","indianred2", "magenta", "green2", "blue2","grey50", "darkviolet","lightblue","chocolate4","navy", "yellow"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K13")

#K=14
temp_dataset_K14<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.14.Q"))
datasetK14<-t(temp_dataset_K14)
colnames(datasetK14)<-angsd_order
datasetK14_reordered<-datasetK14[,paper_order]
#barplot(datasetK14,  col=c("navy", "green2", "cyan1","magenta", "grey50", "blue2", "lightpink","indianred2", "darkviolet","yellow","cornflowerblue","lightblue", "chocolate4", "deepskyblue"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K14")
barplot(datasetK14_reordered,  col=c("navy", "green2", "cyan1","magenta", "grey50", "blue2", "lightpink","indianred2", "darkviolet","yellow","cornflowerblue","lightblue", "chocolate4", "deepskyblue"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K14")

#K=15
temp_dataset_K15<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.15.Q"))
datasetK15<-t(temp_dataset_K15)
colnames(datasetK15)<-angsd_order
datasetK15_reordered<-datasetK15[,paper_order]
#barplot(datasetK15,  col=c("lightpink", "magenta", "deepskyblue","cornflowerblue", "yellow", "darkviolet", "green2","chocolate4", "navy","indianred2","grey50","firebrick2","lightblue", "cyan1", "blue2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K15")
barplot(datasetK15_reordered,  col=c("lightpink", "magenta", "deepskyblue","cornflowerblue", "yellow", "darkviolet", "green2","chocolate4", "navy","indianred2","grey50","firebrick2","lightblue", "cyan1", "blue2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K15")

#K=16
temp_dataset_K16<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.16.Q"))
datasetK16<-t(temp_dataset_K16)
colnames(datasetK16)<-angsd_order
datasetK16_reordered<-datasetK16[,paper_order]
#barplot(datasetK16,  col=c("yellow", "firebrick2", "lightblue","darkviolet", "indianred2", "magenta", "chocolate4","blue2", "lightpink","navy","deepskyblue","green2","black", "cyan1", "grey50", "ivory2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K16")
barplot(datasetK16_reordered,  col=c("yellow", "firebrick2", "lightblue","darkviolet", "indianred2", "magenta", "chocolate4","blue2", "lightpink","navy","deepskyblue","green2","black", "cyan1", "grey50", "ivory2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K16")

#K=17
temp_dataset_K17<-as.matrix(read.table(file = "merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001.admix.17.Q"))
datasetK17<-t(temp_dataset_K17)
colnames(datasetK17)<-angsd_order
datasetK17_reordered<-datasetK17[,paper_order]
#barplot(datasetK17,  col=c("navy", "firebrick2", "yellow","cornflowerblue", "green2", "deepskyblue", "ivory2","blue2", "magenta","indianred2","lightblue","orange","grey50", "darkviolet", "chocolate4", "lightpink", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K17")
barplot(datasetK17_reordered,  col=c("navy", "firebrick2", "yellow","cornflowerblue", "green2", "deepskyblue", "ivory2","blue2", "magenta","indianred2","lightblue","orange","grey50", "darkviolet", "chocolate4", "lightpink", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K17")

dev.off()

#plot the Frobenius error and likelihoods for the different K values
llh<-read.table(file="likelihoods_pcangsd.txt", header = T, dec = ".")
pdf(file = "PCAngsd_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_PCA_Loglikelihood.pdf", width = 11.693, height=8.268)
plot(x=llh$K, y=llh$Frobenius_error, pch=16, cex=1.3, xlab = "K values", ylab = "Frobenius error", col="black", cex.lab=1.2, cex.axis=1.2)
plot(x=llh$K, y=llh$log_likelihood, pch=16, cex=1.3, xlab = "K values", ylab = "Log-likelihood", col="black", cex.lab=1.2, cex.axis=1.2)
dev.off()


##################### NGSAdmix 125inds All individuals (all species) ##############################

#NGSAdmix was run for K values 2 to 17
###Plot the output of the various K values for 125 individuals

#create a pdf file to save the plots
pdf(file = "NGSAdmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_K2_K17_ADMIX_REORDERED.pdf", width = 11.693, height=8.268)

#using the same column names from the pcangsd plots, reorder the individuals in the plots
#read again if necessary (if all variables deleted after finishing the pcangsd plots for example)
#angsd_order<-c(1:125)
#paper_order<-c(c(1:44),c(49:75),c(45:48),c(76,77),c(104:112),c(100:103),c(78:99),c(113:125))


#K=2
ngs_temp_dataset_K2<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_2.qopt"))
ngs_datasetK2<-t(ngs_temp_dataset_K2)
colnames(ngs_datasetK2)<-angsd_order
ngs_datasetK2_reordered<-ngs_datasetK2[,paper_order]
#barplot(ngs_datasetK2,  col=c("grey", "black"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K2")
barplot(ngs_datasetK2_reordered,  col=c("grey", "black"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K2")

#K=3
ngs_temp_dataset_K3<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_3.qopt"))
ngs_datasetK3<-t(ngs_temp_dataset_K3)
colnames(ngs_datasetK3)<-angsd_order
ngs_datasetK3_reordered<-ngs_datasetK3[,paper_order]
#barplot(ngs_datasetK3,  col=c("darkgrey", "black", "ivory2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K3")
barplot(ngs_datasetK3_reordered,  col=c("darkgrey", "black", "ivory2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K3")

#K=4
ngs_temp_dataset_K4<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_4.qopt"))
ngs_datasetK4<-t(ngs_temp_dataset_K4)
colnames(ngs_datasetK4)<-angsd_order
ngs_datasetK4_reordered<-ngs_datasetK4[,paper_order]
#barplot(ngs_datasetK4,  col=c("deepskyblue", "lightgreen", "black", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K4")
barplot(ngs_datasetK4_reordered,  col=c("deepskyblue", "lightgreen", "black", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K4")

#K=5
ngs_temp_dataset_K5<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_5.qopt"))
ngs_datasetK5<-t(ngs_temp_dataset_K5)
colnames(ngs_datasetK5)<-angsd_order
ngs_datasetK5_reordered<-ngs_datasetK5[,paper_order]
#barplot(ngs_datasetK5,  col=c("lightgreen", "black", "deepskyblue", "darkviolet", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K5")
barplot(ngs_datasetK5_reordered,  col=c("lightgreen", "black", "deepskyblue", "darkviolet", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K5")

#K=6
ngs_temp_dataset_K6<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_6.qopt"))
ngs_datasetK6<-t(ngs_temp_dataset_K6)
colnames(ngs_datasetK6)<-angsd_order
ngs_datasetK6_reordered<-ngs_datasetK6[,paper_order]
#barplot(ngs_datasetK6,  col=c("navy", "deepskyblue", "darkviolet", "lightgreen", "black", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K6")
barplot(ngs_datasetK6_reordered,  col=c("navy", "deepskyblue", "darkviolet", "lightgreen", "black", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K6")

#K=7
ngs_temp_dataset_K7<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_7.qopt"))
ngs_datasetK7<-t(ngs_temp_dataset_K7)
colnames(ngs_datasetK7)<-angsd_order
ngs_datasetK7_reordered<-ngs_datasetK7[,paper_order]
#barplot(ngs_datasetK7,  col=c("darkviolet", "black", "deepskyblue", "cyan1", "lightpink", "lightgreen", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K7")
barplot(ngs_datasetK7_reordered,  col=c("darkviolet", "black", "deepskyblue", "cyan1", "lightpink", "lightgreen", "indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K7")

#K=8
ngs_temp_dataset_K8<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_8.qopt"))
ngs_datasetK8<-t(ngs_temp_dataset_K8)
colnames(ngs_datasetK8)<-angsd_order
ngs_datasetK8_reordered<-ngs_datasetK8[,paper_order]
#barplot(ngs_datasetK8,  col=c("lightpink", "darkviolet", "black", "cyan1", "navy", "deepskyblue","lightgreen","indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K8")
barplot(ngs_datasetK8_reordered,  col=c("lightpink", "darkviolet", "black", "cyan1", "navy", "deepskyblue","lightgreen","indianred2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K8")

#K=9
ngs_temp_dataset_K9<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_9.qopt"))
ngs_datasetK9<-t(ngs_temp_dataset_K9)
colnames(ngs_datasetK9)<-angsd_order
ngs_datasetK9_reordered<-ngs_datasetK9[,paper_order]
#barplot(ngs_datasetK9,  col=c("firebrick2","lightgreen", "lightpink", "indianred2", "cyan1", "darkviolet","black","deepskyblue", "blue2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K9")
barplot(ngs_datasetK9_reordered,  col=c("firebrick2","lightgreen", "lightpink", "indianred2", "cyan1", "darkviolet","black","deepskyblue", "blue2"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K9")

#K=10
ngs_temp_dataset_K10<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_10.qopt"))
ngs_datasetK10<-t(ngs_temp_dataset_K10)
colnames(ngs_datasetK10)<-angsd_order
ngs_datasetK10_reordered<-ngs_datasetK10[,paper_order]
#barplot(ngs_datasetK10,  col=c("blue2","lightblue", "indianred2", "magenta", "deepskyblue", "darkviolet","lightpink","firebrick2", "lightgreen", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K10")
barplot(ngs_datasetK10_reordered,  col=c("blue2","lightblue", "indianred2", "magenta", "deepskyblue", "darkviolet","lightpink","firebrick2", "lightgreen", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K10")

#K=11
ngs_temp_dataset_K11<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_11.qopt"))
ngs_datasetK11<-t(ngs_temp_dataset_K11)
colnames(ngs_datasetK11)<-angsd_order
ngs_datasetK11_reordered<-ngs_datasetK11[,paper_order]
#barplot(ngs_datasetK11,  col=c("chocolate4","deepskyblue", "indianred2", "lightblue","magenta", "lightpink", "cyan1", "lightgreen", "darkviolet","blue2", "cornflowerblue"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K11")
barplot(ngs_datasetK11_reordered,  col=c("chocolate4","deepskyblue", "indianred2", "lightblue","magenta", "lightpink", "cyan1", "lightgreen", "darkviolet","blue2", "cornflowerblue"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K11")

#K=12
ngs_temp_dataset_K12<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_12.qopt"))
ngs_datasetK12<-t(ngs_temp_dataset_K12)
colnames(ngs_datasetK12)<-angsd_order
ngs_datasetK12_reordered<-ngs_datasetK12[,paper_order]
#barplot(ngs_datasetK12,  col=c("lightblue","chocolate4", "indianred2", "darkviolet","magenta", "cyan1", "lightpink", "blue2", "cornflowerblue","navy", "deepskyblue", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K12")
barplot(ngs_datasetK12_reordered,  col=c("lightblue","chocolate4", "indianred2", "darkviolet","magenta", "cyan1", "lightpink", "blue2", "cornflowerblue","navy", "deepskyblue", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K12")

#K=13
ngs_temp_dataset_K13<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_13.qopt"))
ngs_datasetK13<-t(ngs_temp_dataset_K13)
colnames(ngs_datasetK13)<-angsd_order
ngs_datasetK13_reordered<-ngs_datasetK13[,paper_order]
#barplot(ngs_datasetK13,  col=c("ivory2","firebrick2","cornflowerblue", "lightpink","deepskyblue", "lightgreen", "navy","darkviolet", "black","lightblue", "chocolate4", "indianred2", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K13")
barplot(ngs_datasetK13_reordered,  col=c("ivory2","firebrick2","cornflowerblue", "lightpink","deepskyblue", "lightgreen", "navy","darkviolet", "black","lightblue", "chocolate4", "indianred2", "cyan1"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K13")

#K=14
ngs_temp_dataset_K14<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_14.qopt"))
ngs_datasetK14<-t(ngs_temp_dataset_K14)
colnames(ngs_datasetK14)<-angsd_order
ngs_datasetK14_reordered<-ngs_datasetK14[,paper_order]
#barplot(ngs_datasetK14,  col=c("lightgreen","navy","firebrick2", "blue3","lightblue","deepskyblue", "indianred1","cornflowerblue", "darkgrey","ivory2", "cyan1", "black", "darkviolet", "lightpink"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K14")
barplot(ngs_datasetK14_reordered,  col=c("lightgreen","navy","firebrick2", "blue3","lightblue","cornflowerblue", "indianred1","deepskyblue", "darkgrey","ivory2", "cyan1", "black", "darkviolet", "lightpink"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K14")

#K=15
ngs_temp_dataset_K15<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_15.qopt"))
ngs_datasetK15<-t(ngs_temp_dataset_K15)
colnames(ngs_datasetK15)<-angsd_order
ngs_datasetK15_reordered<-ngs_datasetK15[,paper_order]
#barplot(ngs_datasetK15,  col=c("navy","ivory2","indianred1", "deepskyblue", "cornflowerblue","lightblue", "firebrick2","chocolate4", "blue2","darkviolet", "magenta", "lightgreen", "cyan1", "lightpink", "orange"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K15")
barplot(ngs_datasetK15_reordered,  col=c("navy","ivory2","indianred1", "deepskyblue", "cornflowerblue","lightblue", "firebrick2","chocolate4", "blue2","darkviolet", "magenta", "lightgreen", "cyan1", "lightpink", "orange"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K15")

#K=16
ngs_temp_dataset_K16<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_16.qopt"))
ngs_datasetK16<-t(ngs_temp_dataset_K16)
colnames(ngs_datasetK16)<-angsd_order
ngs_datasetK16_reordered<-ngs_datasetK16[,paper_order]
#barplot(ngs_datasetK16,  col=c("magenta","lightgreen","yellow", "navy","firebrick2","deepskyblue", "cornflowerblue","black", "ivory2","indianred2", "chocolate4", "lightblue", "blue2", "cyan1", "darkviolet", "lightpink"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K16")
barplot(ngs_datasetK16_reordered,  col=c("magenta","lightgreen","yellow", "navy","firebrick2","deepskyblue", "cornflowerblue","black", "ivory2","indianred2", "chocolate4", "lightblue", "blue2", "cyan1", "darkviolet", "lightpink"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K16")

#K=17
ngs_temp_dataset_K17<-as.matrix(read.table(file = "NGSadmix_merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing_sc_MAF001_17.qopt"))
ngs_datasetK17<-t(ngs_temp_dataset_K17)
colnames(ngs_datasetK17)<-angsd_order
ngs_datasetK17_reordered<-ngs_datasetK17[,paper_order]
#barplot(ngs_datasetK17,  col=c("firebrick2","cyan1","darkviolet", "lightpink","black","chocolate4", "navy","orange", "ivory2","deepskyblue", "indianred2", "lightblue", "blue2", "darkgrey", "cornflowerblue", "magenta", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K17")
barplot(ngs_datasetK17_reordered,  col=c("firebrick2","cyan1","darkviolet", "lightpink","black","chocolate4", "navy","orange", "ivory2","deepskyblue", "indianred2", "lightblue", "blue2", "darkgrey", "cornflowerblue", "magenta", "lightgreen"), xlab="Individual #", ylab="Ancestry", border=NA, main = "K17")

dev.off()
