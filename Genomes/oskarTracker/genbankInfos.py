#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-21-2017
Modified date: 04-13-2017
Description:
Create a csv file where can see if a genome have file which
can help in the several pipelines :
	- protein.fna will help us to found Oskar
	- genomic.bff file can be directly used for Augustus pipeline
	- genomic.gff file can be used for Augustus pipeline if it is combined with genomic.fna file
	- rna.fna presence are given just for information
'''

import os, sys
import pandas as pd
import re
from optparse import OptionParser


def getGenomeInfos(genomeSearch, genomeFolder, folder, currentG, GCA_number, infoDict):
	f = open(genomeSearch,"r")
	oskarGenomes = [x.strip() for x in f.readlines()]
	f.close()

	infoDict[currentG] = {}
	infoDict[currentG]['genome_id'] = GCA_number

	hasToBeChecked = re.findall(r'(^[A-Z]{3}_[0-9]{9}).[0-9]{1}', currentG)[0]

	for genome in oskarGenomes :
		checked = re.findall(r'(^[A-Z]{3}_[0-9]{9}).[0-9]{1}', genome)[0]
		if hasToBeChecked == checked :
			print("Oskar already identified in "+ hasToBeChecked)
			infoDict[currentG]['oskar'] = 'Yes'

	if 'Oskar presence' not in infoDict[currentG]:
		infoDict[currentG]['oskar'] = 'N/A'

def getallGenomesInfos(pathcsv, genbankFolder) :
	df = pd.read_csv(pathcsv)
	df_i = df.sort_values(["order_id", "family_id"], ascending=[0,1])[['order_name', 'family_name', 'sp_name', 'genome_id']]

	allGenomesInfos = []

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

			hasToBeFound = re.findall(r'(^[A-Z]{3}_[0-9]{9}.[0-9]{1})', current)[0]

			for genome in df_i['genome_id'] :
				regExpr = re.findall(r'^([A-Z]{3}_[0-9]{9}.[0-9]{1})', genome)[0]
				if hasToBeFound == regExpr :
					allGenomesInfos.append(getGenomeInfos(genomeSearch, genomeFolder, folder, current, hasToBeFound, {}))

	return df_i, allGenomesInfos

def buildFile(allGenomesInfos, df_i):
	inf_df = []

	for elem in allGenomesInfos :
		for key in elem.keys():
			inf_df.append(elem[key])

	df_m = pd.DataFrame(inf_df)
	df_m = df_m[['genome_id','oskar']]

	df_f = pd.merge(df_i, df_m, how='inner', on=['genome_id'], sort=False, suffixes=('_x', '_y'), copy=True, indicator=False)
	df_f = df_f.set_index('sp_name')

	df_f.to_csv('insect_oskarInfos.csv', sep=', ', na_rep='', header=True)


# Main

parser = OptionParser()
parser.add_option("-g", "--genbank_folder", dest="genbankFolder", default="None",
                  help="[Required] Location of the folder containing Genbank genomes")
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file curated_InsectGenomesInfos.csv")
parser.add_option("-s", "--skipped_genomes", dest="genomeSearch", default="None",
                  help="[Required] Location of the file which contains all the genomes where we already know that Oskar gene exist")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv
genbankFolder = options.genbankFolder
refseqFolder = options.refseqFolder
genomeSearch = options.genomeSearch

if genbankFolder == "None":
	print("Genbank Folder must be provided")
	sys.exit(1)

if pathcsv == "None":
	print("curated_InsectGenomesInfos.csv file must be provided")
	sys.exit(1)

if genomeSearch == "None":
	print("GenomeSearch file must be provided")
	sys.exit(1)

pre_insectDataFrame, allGenomesInfos = getAllGenomesInfos(pathcsv, genbankFolder, refseqFolder)
print("Traking information terminated")
buildFile(allGenomesInfos, pre_insectDataFrame)
