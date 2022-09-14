#$ -S /bin/bash
#$ -N psi4_dft
#$ -R y
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_rss=24G
#$ -q cpus
#$ -l cpus=1
#$ -pe smp 8
# -fin *.py
# -fin b3lyp
# -fin *.xyz
# -fin *.molden
# -fin *.json

source /home/crduan/.bashrc
conda activate /home/crduan/miniconda/envs/mols_py36
export PSI_SCRATCH='./'

echo 'psi4 scr: ' $PSI_SCRATCH
python -u loop_run.py  > $SGE_O_WORKDIR/nohup1.out
python -u loop_run.py  > $SGE_O_WORKDIR/nohup2.out
python -u loop_run.py  > $SGE_O_WORKDIR/nohup3.out

mv * $SGE_O_WORKDIR
