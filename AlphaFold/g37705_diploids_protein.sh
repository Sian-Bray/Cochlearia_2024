#!/bin/bash
#PBS -N g37705_diploids_protein
#PBS -q gpu@meta-pbs.metacentrum.cz
#PBS -l select=1:ncpus=8:ngpus=1:mem=200gb:scratch_local=10gb:gpu_cap=cuda75
#PBS -l walltime=24:00:00
#PBS -m ae

# !!! gpu_cap=cuda80 == only ZIA cluster, see https://wiki.metacentrum.cz/wiki/GPU_clusters


#############################################
# !!change these paths according your dirs !!!
#############################################
INPUTS_DIR=/storage/praha1/home/razi/Sian_thank_you_Razi  #this DIR is for multiple inpuut files, see lines 20-21
OUTPUTS_DIR=/storage/praha1/home/razi/Sian_thank_you_Razi/AlphaFold_outputs/outputs-$PBS_JOBID
#############################################

mkdir -p $SCRATCHDIR/tmp  # tmp directory, some jobs needs >1GB

#FASTA_INPUTS=`ls -1 $INPUTS_DIR/*.fasta | xargs -n 1 basename | sed -e "s/^/\/fasta.inputs\//g" -e "s/^$//g"| tr "\n" "," | sed -e "s/,$//g"`  # this line takes all *.fasta files in $INPUT_DIR
FASTA_INPUTS=/storage/praha1/home/razi/Sian_thank_you_Razi/AlphaFold_fastas/All_Proteins/g37705_diploids_protein.fasta
echo --$FASTA_INPUTS--

HOSTNAME=`hostname`
if [[ $HOSTNAME =~ adan.* ]]; then
    STORAGE=/storage/vestec1-elixir/projects
else
    STORAGE=/storage/brno11-elixir/projects
fi

export TMPDIR=/scratch
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=10.0

singularity exec --nv -B $STORAGE/alphafold/alphafold.db/:/alphafold.db --pwd /app/alphafold -B $SCRATCHDIR:/scratch -B /$SCRATCHDIR/tmp:/tmp -B $INPUTS_DIR:/fasta.inputs\
  /storage/brno11-elixir/projects/alphafold/AF2.1-ubuntu20.04-cuda11.1.sif  \
      python /app/alphafold/run_alphafold.py \
        --fasta_paths=$FASTA_INPUTS \
        --is_prokaryote_list=false \
        --max_template_date=2020-05-14 \
        --model_preset=monomer \
        --data_dir=/alphafold.db \
        --output_dir=/scratch  \
        --bfd_database_path=/alphafold.db/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt  \
        --pdb70_database_path=/alphafold.db/pdb70/pdb70 \
        --uniref90_database_path=/alphafold.db/uniref90/uniref90.fasta \
        --uniclust30_database_path=/alphafold.db/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
        --mgnify_database_path=/alphafold.db/mgnify/mgy_clusters.fa \
        --template_mmcif_dir=/alphafold.db/pdb_mmcif.2.1/mmcif_files \
        --obsolete_pdbs_path=/alphafold.db/pdb_mmcif.2.1/obsolete.dat \
        --db_preset=full_dbs

mkdir -p $OUTPUTS_DIR
cp -r  $SCRATCHDIR/* $OUTPUTS_DIR/ &&  rm -r $SCRATCHDIR/*
