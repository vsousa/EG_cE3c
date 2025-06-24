#!/bin/bash

#S. L. Mendes 20/08/2021
#Based on D. Frei script (https://github.com/freidavid/Genomic-Consequences-of-Speciation-Reversal/blob/main/lsf_submit_scripts/merging_single_beagle_min30.lsf)

#SETTINGS
#Name of the four beagle files that we want to merge (common part of the name)
NAME="lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing"
#Tag of the reference genome used
REFTAG="sc"
#Name of the final beagle file:
MERGED="merged_lcwgs_squalius_125inds_minIndDP2_maxIndDP14_1pcmissing"
#Path to ANGSD output folder, that is, where the four output beagle files are and 
#where we also want to save the merged beagle files
ANGSDFOLDER=/path/to/output/angsd/folder;
#Path to text file with the distinct part of the names of the beagle files we want to merge
TERM_LIST=/media/shared/smendes/8_GL_minIndDP2_maxIndDP2SD_125inds_scephalus/chromossome_names_scephalus.txt;


#Header taken from the first (_1) angsd output beagle file that will be merged later. Files 
#to be merged will be attached to this header.
zcat ${ANGSDFOLDER}/${NAME}_CM040792.1:_${REFTAG}.beagle.gz | head -n1 > ${ANGSDFOLDER}/${MERGED}_${REFTAG}.beagle

#Attach all the beagle files to the new merged one
while read line 
do
	zcat ${ANGSDFOLDER}/${NAME}_${line}_${REFTAG}.beagle.gz | sed -n '1!p' >> ${ANGSDFOLDER}/${MERGED}_${REFTAG}.beagle
done < ${TERM_LIST}

#gzip the merged beagle file
gzip ${ANGSDFOLDER}/${MERGED}_${REFTAG}.beagle

