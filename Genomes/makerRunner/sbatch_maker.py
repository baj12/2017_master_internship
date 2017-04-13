#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04/03/2017
Modified date: 04/04/2017

Description:
Create sbatch files for each representant order
which be after running to execute Maker
'''


import os, sys
import re

for folder in os.listdir('/n/extavour_lab/Lab/Everyone/savy/Maker'):
        f = open('{}/{}/{}.sh'.format('/n/extavour_lab/Lab/Everyone/savy/Maker',folder,folder),'w')
        f.write('#!/bin/bash\n#\n')
        f.write('#SBATCH -J {} \n'.format(folder))
        f.write('''#SBATCH -n 32 # Number of cores
#SBATCH -p general # Partition
#SBATCH --mem 16000 # Memory request
#SBATCH -o maker.out # Standard output
#SBATCH -e maker.err # Standard error
#
#
#
module load gcc/5.2.0-fasrc01 openmpi/2.0.1-fasrc01 maker/2.31.8-fasrc01
module load RepeatMasker/4.0.5-fasrc04
module load exonerate/2.4.0-fasrc01
''')
    	f.write('maker --fix_nucleotides')
        f.write('\n')
        f.close()
