#!/bin/bash
source ~/venv/optim2updt/bin/activate

cd /home/s/slautenbach/projects/optim/LPJ_pareto/code/version_fpu_updated_dec_2019/post_processing
python2.7 createTableFromFront_main_fpu.py 
#python2.7 createGenomeTabelesFPU_main.py