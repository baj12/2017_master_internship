#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04/17/2017
Modified date: 04/04/2017

Description:
Create sbatch files for each genome folder
which be after running to execute augustus training
'''


import os, sys, re
from Bio import Entrez
Entrez.email = 'savandara.besse@gmail.com'


path = '/n/regal/extavour_lab/savy/augustus_training'
config_path = '/n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/config/species'
for folder in os.listdir(path):
    genome_id = re.findall(r'(^[A-Z]{3}_[0-9]{9}.[0-9]{1})',folder)[0]
    search_handle = Entrez.esearch(db="genome", term=genome_id)
    orgn_id = Entrez.read(search_handle)['IdList']

    handle = Entrez.esummary(db="genome", id=orgn_id)
    name = Entrez.read(handle)[0]['Organism_Name'].replace(' ', '_')

    f = open('{}/{}_training.sh'.format(path,name),'w')
    f.write('#!/bin/bash\n#\n')
    f.write('#SBATCH -J {}_training \n#SBATCH -o {}_training.out # Standard output\n#SBATCH -e {}_training.err # Standard error\n'.format(name,name,name))
    f.write('''#SBATCH -n 1 # Number of cores
#SBATCH -p general # Partition
#SBATCH --mem 2000 # Memory request
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#
#
#
module load augustus/3.0.3-fasrc02
module load boost/1.55.0-fasrc01
module load bamtools/2.3.0-fasrc01
module load samtools/0.1.19-fasrc01
module load bcftools/1.0-fasrc01
module load htslib/1.1-fasrc01
module load zlib/1.2.8-fasrc02
module load tabix/0.2.6-fasrc01
''')
    f.write('/n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/scripts/gff2gbSmallDNA.pl {}/{}_genomic.gff {}/{}_genomic.fna 5000 {}/{}_genes.gb \n'.format(path,folder,path,folder,path,folder))
    f.write('/n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/scripts/randomSplit.pl {}/{}_genes.gb 100\n'.format(path,folder))
    f.write('/n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/scripts/randomSplit.pl {}/{}_genes.gb 100\n'.format(path,folder))
    f.write('/n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/scripts/new_species.pl --species={} --AUGUSTUS_CONFIG_PATH={}\n'.format(name,config_path))
    f.write('/n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/bin/etraining --species={} --AUGUSTUS_CONFIG_PATH={} {}/{}_genes.gb.train\n'.format(name,config_path,path,folder))
    f.write('/n/sw/fasrcsw/apps/Core/augustus/3.0.3-fasrc02/bin/augustus --species={} --AUGUSTUS_CONFIG_PATH={} {}/{}_genes.gb.test\n'.format(name,config_path,path,folder))
    f.close()
