#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-26-2017
Modified date: 03-27-2017

Description:
Create a csv file where can see for each organism and tell us :
- its order and its family
- if they have protein annotations file, its protein number
- if they have gene annotations file, its gene number
- if they have est in NCBI, its EST number

'''

import os, sys
import pandas as pd
import re
from optparse import OptionParser

def openGenomeFolder(genomeFolder, genomeID):
    genomeidFolder = os.listdir(genomeFolder)
	for folder in genomeidFolder :
		path = os.path.join(genomeFolder, folder)
		if not os.path.isfile(path) :
			genomeVersion = os.listdir(path)

            for genome in genomeVersion:
                if genomeID in genome :
                    path = os.path.join(genomeFolder, folder, genome)
                    print(path)


def countProtein(pathcsv):
    df = pd.read_csv(pathcsv)

    for species in df.index :
        if df["Protein annotations"][species] == "Yes" :
            current_genomeID = df["Accession number"]
            openGenomeFolder(genomeFolder, current_genomeID)



parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file curated_InsectGenomesInfos.csv")


(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
pathcsv = options.pathcsv

if genomeFolder == "None":
	print("Genome Folder must be provided")
	sys.exit(1)

if pathcsv == "None":
	print("curated_InsectGenomesInfos.csv file must be provided")
	sys.exit(1)

countProtein(csv)
