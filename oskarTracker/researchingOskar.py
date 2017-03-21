#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-17-2017
Modified date: 03-20-2017 
Description: 

Skipped genomes which are already know to have the Oskar sequence
If it is not the case, the script write the command to execute Exonerate by itself
At the end, you will obtain exonerate output, which will tell you if you have an oskar sequence or not in your genome
'''

import os, sys
import re
from optparse import OptionParser
from subprocess import Popen
import subprocess
import time


def skipped_Genomes(genomeSearch):
	f = open(genomeSearch,"r")
	oskarGenomes = [x.strip() for x in f.readlines()]
	f.close()

	return oskarGenomes


def add_Command(genomeFolder, folder, currentG, exonerateCMD):
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
		exonerateCMD = ["/mirror/bin/exonerate-2.2.0/bin/exonerate", "--showtargetgff", "yes",  "--fsmmemory", "2500", "--model", "protein2genome",  "--percent", "50", "--target", path,  "--query", "OSKAR_FINAL.fasta"] 

	return exonerateCMD

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
					print(i, current + " already known to have Oskar sequence")
				else :
					cmdToRun.append(add_Command(genomeFolder, folder, current, []))

	return genomeList, cmdToRun

def run_AllCmd(cmdToRun) :
	n = 0
	running = {}
	process = {}
	name = {}
	path = {}
	UID = 0
	log = open('oskar_tracker.log','w')

	while cmdToRun:
		try:
			while n < nb_CPU:
				UID += 1
				newcmd = cmdToRun.pop()
				print(" ".join(newcmd))
				P = Popen(newcmd)
				print("Launched {}/{}".format(UID, len(cmdToRun)))
				running[UID] = 1
				process[UID] = P
				name[UID] = " ".join(newcmd)
				n += 1
		
			for UID in running.keys():
				if process[UID].poll() == 0:
					del running[UID]
					print("FINISHED" + name[UID])
					n -= 1
			time.sleep(1)
		except KeyboardInterrupt:
			for k in running.keys():
				os.remove(path[UID])
			print("cleaned running file")
			print("exiting")


def remove_Genome(genomeFolder, genomeList):
	print("Removing Genome has started")

	genomeidFolder = os.listdir(genomeFolder)
	for folder in genomeidFolder :
		for genome in genomeList : 
			path = os.path.join(genomeFolder, folder, genome, genome+"_genomic.fna")
			if os.path.isfile(path):
				os.system("rm -rf "+path)

	print("AllGenomes are removed")


# Main

parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-p", "--file", dest="genomeSearch", default="None",
                  help="[Required] Location of the file which contains all the genomes where we already know that Oskar gene exist")
parser.add_option("-c", "--cpu", dest="cpu", default=1,
                  help="Number of CPU to use for the analysis. Each core is used to launch one instance of SNAP or Augustus.")

(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
genomeSearch = options.genomeSearch

if genomeFolder == "None":
	print("Genome Folder must be provided")
	sys.exit(1)

if genomeSearch == "None":
	print("GenomeSearch file must be provided")
	sys.exit(1)

try:
	nb_CPU = int(options.cpu)
except:
	print("Wrong number of CPUs !")
	sys.exit(1)

oskarGenomes = skipped_Genomes(genomeSearch)
genomeList, cmdToRun = build_cmdToRun(genomeFolder, oskarGenomes)
run_AllCmd(cmdToRun)
remove_Genome(genomeFolder, genomeList)