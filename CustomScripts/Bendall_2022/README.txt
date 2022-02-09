Scripts from Bendall et al (2022)

This contains the following folders:
- Scripts_run_SLiM: input files and scripts to run SLIM for simulation study under IM model with symmetric migration
- Scripts_process_SLiM: scripts to process SLiM output files for simulation study under IM model with symmetric migration
- Scripts_SLiM_Neodiprion: input files and scripts to run and process SLIM simulations under demographic history inferred for Neodiprion sawflies
- Scripts_figures_SLiM_summaryFiles: scrtips to produce the figures in the manuscript. For figures about sawflies check Scripts_SLiM_Neodiprion/Scripts_plots_FiguresManuscript

SIMULATION STUDY
Simulations were done using SLiM version 3.3 (built Aug  5 2019 18:29:03).
Check more details in README files in folders:
 ./Scripts_run_SLiM
 ./Scripts_run_SLiM/scaled
Within those folders you can find the SLiM template input files, ending with *.s and the bash scripts used to modify the template files and launch SLiM.

Launch SLIM simulations
Check more details in README file in folder
./Scripts_run_SLiM/scaled

# 1. Create ancestral population define combination of parameters in CreateAncestralPops.sh
./CreateAncestralPops.sh

# 2. Simulate the divergence of populations
#  specify the template file (either sim_divsel3_scaled_new_template.s
#  or sim_divsel3_scaled_template.s)
#  use the same combination of parameters as in CreateAncestralPops.sh
./SimFromAncestralPops.sh

Process output files from SLiM
Check more details in README files in folders: 
./Scripts_process_SLiM_output
./Scripts_process_SLiM_output/scaled

1.	compute statistics for linked selection based on samples of 20 genomes of 500kb from each pop average across the chromosome and window-based analysis read_ms_createSummaryFiles.r
2.	compute statistics for selected site from simulations of chromosomes with 500kb read_traj_createSummaryFiles.r
3.	further process statistics of ms to obtain means and quantiles of pi, dxy and fst process the summary files of ms to create simple tables that can be used for plots. process_ms_summaryFiles.r
4.	further process statistics of ms to obtain means and quantiles of pi and fst along the chromosome  on windows base.  process_ms_summaryFiles_windowscan.r
5.	read ms files to compute LD r2 statistic read_ms_compute_LD.r
6.	Process the LD statistics files process_ld.r


SIMULATIONS ACCORDING TO Neodiprion INFERRED DEMOGRAPHIC HISTORY
Simulations were done using SLiM version 3.3 (built Aug  5 2019 18:29:03).

Check more details in README file in folder:
 ./Scripts_SLiM_Neodiprion
Launch SLIM simulations according to Neodiprion demographic history
Check README file in:
./Scripts_SLiM_Neodiprion/Scripts_Run_SLiM

# 1. Create ancestral population define combination of parameters in CreateAncestralPops.sh
./CreateAncestralPops.sh

# 2. Simulate the divergence of populations 
# specify the template file (either sim_divsel3_xaut_newmut_template.s or 
# sim_divsel3_xaut_template.s) use the same combination of parameters as in
# CreateAncestralPops.sh
./SimFromAncestralPops.sh

Process output files from SLiM according to Neodiprion demographic history
Check more details in README file in folder: 
./Scripts_SLiM_Neodiprion/Scripts_Process_SLiM
1.	compute statistics for linked selection based on samples of 20 genomes from each pop readms_createSummaryFiles_sawflies.r
2.	process statistics of ms to obtain means and quantiles that can be used for plots. process_ms_summaryFiles_neutralSFS_windows.r


FIGURES FOR MANUSCRIPT
Please check the README files and scripts in the following folder:
For simulation study:
./ Scripts_figures_SLiM_summaryFiles

For results related with simulations of Neodiprion sawflies:
./Scripts_SLiM_Neodiprion/Scripts_plots_FiguresManuscript
