#S.L.Mendes
#01/02/2022
#script to calculate abba/baba with angsd per groups
#edited on 10.04.2023 for the Squalius cephalus genome in 5000000 bp windows

#Folder where we will run the D-statistics and where this script MUST be run
DSTATFOLDER=/path/to/dstat/folder
cd ${DSTATFOLDER}

#We need a file with the list of bam files of the individuals we want to test. Note that the last 2 individuals are the outgroup
#Because we want to perform the test using different outgroups (S. torgalensis and S. aradensis), we need two sets of bam files
BAMFILESOUTGT=/media/shared/smendes/10_Dstat_scephalus/New2023data_125inds/18_Lizandro_Pego_Ardila/listBamFiles_outgT.txt
BAMFILESOUTGA=/media/shared/smendes/10_Dstat_scephalus/New2023data_125inds/18_Lizandro_Pego_Ardila/listBamFiles_outgA.txt
#file with the name of each population
#Because we want to perform the test using different outgroups (S. torgalensis and S. aradensis), we need two files
NAMEFILEOUTGT=/media/shared/smendes/10_Dstat_scephalus/New2023data_125inds/18_Lizandro_Pego_Ardila/popNames_outgT.name
NAMEFILEOUTGA=/media/shared/smendes/10_Dstat_scephalus/New2023data_125inds/18_Lizandro_Pego_Ardila/popNames_outgA.name


#file with size of each population
#because the size of our outgroup populations is the same (2 individuals each), we only need one file 
SIZEFILE=/media/shared/smendes/10_Dstat_scephalus/New2023data_125inds/18_Lizandro_Pego_Ardila/sizeFile.size

#file with the 25 well assmembled chormossomes of the S. cephalus genome
REGIONS=/media/shared/smendes/10_Dstat_scephalus/New2023data_125inds/18_Lizandro_Pego_Ardila/chromossome_names_scephalus.txt

#Name/Tag of this test (this identifies the populations used to perform the test in the file names)
TAGNAME=Lizandro_Pego_Ardila

#Tag/Code of the outgroups (this identifies the outgroup used in the file names)
#S. torgalensis as the outgroup
OUTGT=outT
#S. aradensis as the outgroup
OUTGA=outA

#Error file for Rscript 
echo "NA" > ${DSTATFOLDER}/errorFile.error 
echo "NA" >> ${DSTATFOLDER}/errorFile.error 
echo "NA" >> ${DSTATFOLDER}/errorFile.error 
echo "NA" >> ${DSTATFOLDER}/errorFile.error 



###################### S. torgalensis as the outgroup #######################

#Calculate the abba/baba (D-statistic)
angsd -doAbbababa2 1 -out ${DSTATFOLDER}/dstat_${TAGNAME}_${OUTGT} -nThreads 1 \
-doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-setMinDepthInd 2 -setMaxDepthInd 14 \
-bam ${BAMFILESOUTGT} -sizeFile ${SIZEFILE} -rf ${REGIONS} -useLast 1 -blockSize 5000000

#get the readable results using the Rscript
Rscript ./estAvgError.r angsdFile=${DSTATFOLDER}/dstat_${TAGNAME}_${OUTGT} out=${DSTATFOLDER}/dstat_${TAGNAME}_${OUTGT}_final sizeFile=${SIZEFILE} nameFile=${NAMEFILEOUTGT} errFile=${DSTATFOLDER}/errorFile.error > dstat_${TAGNAME}_${OUTGT}_final.log




###################### S. aradensis as the outgroup #######################

#Calculate the abba/baba (D-statistic)
angsd -doAbbababa2 1 -out ${DSTATFOLDER}/dstat_${TAGNAME}_${OUTGA} -nThreads 1 \
-doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-setMinDepthInd 2 -setMaxDepthInd 14 \
-bam ${BAMFILESOUTGA} -sizeFile ${SIZEFILE} -rf ${REGIONS} -useLast 1 -blockSize 5000000

#get the readable results using the Rscript
Rscript ./estAvgError.r angsdFile=${DSTATFOLDER}/dstat_${TAGNAME}_${OUTGA} out=${DSTATFOLDER}/dstat_${TAGNAME}_${OUTGA}_final sizeFile=${SIZEFILE} nameFile=${NAMEFILEOUTGA} errFile=${DSTATFOLDER}/errorFile.error > dstat_${TAGNAME}_${OUTGA}_final.log




