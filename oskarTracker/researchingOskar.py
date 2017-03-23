#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-17-2017
Modified date: 03-21-2017

Description: 
Skipped genomes which are already know to have the Oskar sequence
If it is not the case, the script write the command to execute Exonerate in in a sh_file to be executed by Slurm
At the end, you will obtain exonerate output, which will tell you if you have an oskar sequence or not in your genome
'''

import os, sys
import re
from optparse import OptionParser
from subprocess import Popen
import subprocess


def skipped_Genomes(genomeSearch):
	f = open(genomeSearch,"r")
	oskarGenomes = [x.strip() for x in f.readlines()]
	f.close()

	return oskarGenomes


def add_Command(genomeFolder, folder, currentG, exonerateARG):
	print("Opening Genome: {}".format(currentG))

	path = os.path.join(genomeFolder, folder, currentG, currentG+"_genomic.fna.gz")
	if os.path.isfile(path):
		if not os.path.isfile(os.path.join(genomeFolder, folder, currentG, currentG+"_genomic.fna")):
			P = Popen(['gunzip',path])
			ret = P.wait()
			if ret != 0:
				print("Error Gunzipping !")
	path = os.path.join(genomeFolder, folder, currentG, currentG+"_genomic.fna")
	if os.path.isfile(path):
		exonerateARG[currentG] = {"$1":path} 

	return exonerateARG


def build_cmdToRun(genomeFolder, oskarGenomes):
	cmdToRun = []
	genomeList = []
	genomeidFolder = os.listdir(genomeFolder)
	for folder in genomeidFolder :
		path = os.path.join(genomeFolder, folder)
		if not os.path.isfile(path) :
			genomeVersion = os.listdir(path)

			if len(genomeVersion) == 1 :
				current = genomeVersion[0]

			if len(genomeVersion) > 1 :
				maxVersion = 0
				for genome in genomeVersion:
					nb = re.findall(r'^[A-Z]{3}_[0-9]{9}.([0-9]{1})', genome)[0]
					if int(nb) > maxVersion :
						maxVersion = int(nb)
						lastVersion = genomeVersion.index(genome)
				current = genomeVersion[lastVersion]

			genomeList.append(current)

			hasToBeFound = re.findall(r'(^[A-Z]{3}_[0-9]{9}).[0-9]{1}', current)[0]

			for genome in oskarGenomes :
				found = re.findall(r'(^[A-Z]{3}_[0-9]{9}).[0-9]{1}', genome)[0]
				if hasToBeFound == found :
					print(hasToBeFound+ " already known to have Oskar sequence")
					current = None

			if current is not None :
				cmdToRun.append(add_Command(genomeFolder, folder, current, {}))

	return genomeList, cmdToRun


def run_AllCmd(cmdToRun, proteinDatabase) :
	if not os.path.isdir('SLURM'):
		os.mkdir('SLURM')

	jobID = 1
	for cmd in cmdToRun:
		for genome in cmd.keys() : 
			f = open('./SLURM/exonerate_{}.sh'.format(jobID),'w')
			f.write('#!/bin/bash\n#SBATCH \n\n')
			f.write('\n/mirror/bin/exonerate-2.2.0/bin/exonerate --query ' +proteinDatabase+ ' --target ' +cmd[genome]['$1']+ ' --model protein2genome --percent 50 --showtargetgff yes --showalignment no --showvulgar no -M 1500 \n')
			f.close()
		jobID += 1

	f = open('SLURM_exonerate.sh','w')
	f.write('#!/bin/bash\n#SBATCH -J tophat\n#SBATCH --mem 2800\n#SBATCH -n 1\n#SBATCH -e log_exonerate.err\n\n')
	f.write('\n\nbash {}/SLURM/exonerate_"${{SLURM_ARRAY_TASK_ID}}".sh'.format(os.path.abspath('.')))
	f.write('\n\n# To run please execute\n# sbatch --array=1-{} SLURM_exonerate.sh'.format(len(cmdToRun)))
	f.close()


# Main

parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-s", "--skipped_genomes", dest="genomeSearch", default="None",
                  help="[Required] Location of the file which contains all the genomes where we already know that Oskar gene exist")
parser.add_option("-d","--database", dest="proteinDatabase", default="None",
                  help="[Required] Location of the file which contains uniprot_AllProteins and Oskar sequences")

(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
genomeSearch = options.genomeSearch
proteinDatabase = options.proteinDatabase

if genomeFolder == "None":
	print("Genome Folder must be provided")
	sys.exit(1)

if genomeSearch == "None":
	print("GenomeSearch file must be provided")
	sys.exit(1)

if proteinDatabase == "None":
	print("proteinDatabase file must be provided")
	sys.exit(1)

oskarGenomes = skipped_Genomes(genomeSearch)
genomeList, cmdToRun = build_cmdToRun(genomeFolder, oskarGenomes)
run_AllCmd(cmdToRun, proteinDatabase)