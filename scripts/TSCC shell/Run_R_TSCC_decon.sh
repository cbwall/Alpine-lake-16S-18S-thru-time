#!/bin/bash
#PBS -q hotel
#PBS -N 18S_dada_decon
#PBS -l nodes=1:ppn=2
#PBS -l walltime=01:00:00
#PBS -o /projects/ps-shurinlab/users/cbwall/Yos_water_18S/output/outs/Rscript_out
#PBS -e /projects/ps-shurinlab/users/cbwall/Yos_water_18S/output/errs/Rscript_err
#PBS -V
#PBS -M cbwall@ucsd.edu
#PBS -m ae
#PBS -A shurin-group

export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH
module load R/4.1.2
Rscript --no-save /projects/ps-shurinlab/users/cbwall/scripts/18S_dada_decon.R

