#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-21-2017
Modified date: 03-21-2017 
Description: 
Create a csv file where can see if a genome have or not annotations
'''

import os, sys
import re
from optparse import OptionParser

def getInfosGenome(genomeFolder, folder, currentG):
	print("Opening Genome: {}".format(currentG))
	path = os.path.join(genomeFolder, folder, currentG)
	allInfos = os.listdir(path)

	return allInfos


def openFolder(genomeFolder) :
	genomeidFolder = os.listdir(genomeFolder)
	for folder in [genomeidFolder[0]] :
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

	return getInfosGenome(genomeFolder, folder, current)


# Main

parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")

(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder

result = openFolder(genomeFolder)
print(result)
