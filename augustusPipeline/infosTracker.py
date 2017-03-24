#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-21-2017
Modified date: 03-23-2017 
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


def makeDecisions(infoDict, currentG):

	print(currentG+ " tracking information in progress")

	if infoDict[currentG]["Oskar presence"] == 'N/A' :
		if infoDict[currentG]["Protein annotations"] == 'Yes' :
			infoDict[currentG]["TO DO"] = 'Find Oskar'
		elif infoDict[currentG]["GeneBank annotations"] == 'Yes':
			infoDict[currentG]["TO DO"] = 'Check GFF3 annotations'
		elif infoDict[currentG]["GeneBank annotations"] == 'No' and (infoDict[currentG]["GFF3 annotations"] == 'Yes' and infoDict[currentG]["Genomic fasta file"] == 'Yes'):
			infoDict[currentG]["TO DO"] = 'Check GFF3 annotations'
		elif infoDict[currentG]["GFF3 annotations"] == 'No' and infoDict[currentG]["Genomic fasta file"] == 'Yes':
			infoDict[currentG]["TO DO"] = 'Use GeneMark + Exonerate + Evidence Modeller + Augustus'
		elif infoDict[currentG]["Genomic fasta file"] == 'No' :
			infoDict[currentG]["TO DO"] = 'Should look to the folder '+ currentG
		else : 
			infoDict[currentG]["TO DO"] = "Good question ..."

	elif infoDict[currentG]["Oskar presence"] == 'Yes' :
		infoDict[currentG]["TO DO"] = 'Oskar already known in this genome'

	return infoDict


def getGenomeInfos(genomeSearch, genomeFolder, folder, currentG, GCA_number, infoDict):
	f = open(genomeSearch,"r")
	oskarGenomes = [x.strip() for x in f.readlines()]
	f.close()

	infoDict[currentG] = {}
	infoDict[currentG]['Accession number'] = GCA_number

	hasToBeChecked = re.findall(r'(^[A-Z]{3}_[0-9]{9}).[0-9]{1}', currentG)[0]

	for genome in oskarGenomes :
		checked = re.findall(r'(^[A-Z]{3}_[0-9]{9}).[0-9]{1}', genome)[0]
		if hasToBeChecked == checked :
			print("Oskar already identified in "+ hasToBeChecked)
			infoDict[currentG]['Oskar presence'] = 'Yes'
	
	if 'Oskar presence' not in infoDict[currentG]:
		infoDict[currentG]['Oskar presence'] = 'N/A'
			
	path = os.path.join(genomeFolder, folder, currentG)
	allInfos = os.listdir(path)

	for info in allInfos :
		if '_protein.faa.gz' in info :
			infoDict[currentG]['Protein annotations'] = 'Yes'
		if '_genomic.gbff.gz' in info :
			infoDict[currentG]['GeneBank annotations'] = 'Yes'
		if '_genomic.gff.gz' in info :
			infoDict[currentG]['GFF3 annotations'] = 'Yes'
		if '_genomic.fna' in info :
			infoDict[currentG]['Genomic fasta file'] = 'Yes'
		if '_rna.fna.gz' in info :
			infoDict[currentG]['RNA fasta file'] = 'Yes'

	if 'Protein annotations' not in infoDict[currentG]:
		infoDict[currentG]['Protein annotations'] = 'No'
	if 'GeneBank annotations' not in infoDict[currentG]:
		infoDict[currentG]['GeneBank annotations'] = 'No'
	if 'GFF3 annotations' not in infoDict[currentG]:
		infoDict[currentG]['GFF3 annotations'] = 'No'
	if 'Genomic fasta file' not in infoDict[currentG]:
		infoDict[currentG]['Genomic fasta file'] = 'No'
	if 'RNA fasta file' not in infoDict[currentG]:
		infoDict[currentG]['RNA fasta file'] = 'No'

	return makeDecisions(infoDict, currentG)


def getAllGenomesInfos(pathcsv, genomeFolder) :
	df = pd.read_csv(pathcsv)
	df_i = df.sort_values(["Order id", "Family id"], ascending=[0,1])[['Order name', 'Family name', 'Species name', 'Accession number']]

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

			for genome in df_i['Accession number'] :
				regExpr = re.findall(r'^([A-Z]{3}_[0-9]{9}.[0-9]{1})', genome)[0]
				if hasToBeFound == regExpr :
					allGenomesInfos.append(getGenomeInfos(genomeSearch, genomeFolder, folder, current, hasToBeFound, {}))

	print("Traking information terminated")

	return df_i, allGenomesInfos


def buildpipelineGuide(allGenomesInfos, df_i):
	inf_df = []

	for elem in allGenomesInfos : 
		for key in elem.keys():
			inf_df.append(elem[key])

	df_m = pd.DataFrame(inf_df)
	df_m = df_m[['Accession number','Oskar presence','Protein annotations','GeneBank annotations','GFF3 annotations','Genomic fasta file','RNA fasta file','TO DO']]

	df_f = pd.merge(df_i, df_m, how='inner', on=['Accession number'], sort=False, suffixes=('_x', '_y'), copy=True, indicator=False)
	df_f = df_f.set_index('Species name')
	
	f = open('pipelineGuide.csv','w', encoding="utf-8")

	f.write("Species name,")
	for column in df_f.columns:
		if column == "TO DO":
			f.write(column+'\n')
		else :
			f.write(column+',')

	for species in df_f.index:
		f.write(species+',')
		for column in df_f.columns:
			if column == 'TO DO':
				f.write(str(df_f[column][species])+'\n')
			else:
				f.write(str(df_f[column][species])+',')

	f.close()


# Main

parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file curated_InsectGenomesInfos.csv")
parser.add_option("-s", "--skipped_genomes", dest="genomeSearch", default="None",
                  help="[Required] Location of the file which contains all the genomes where we already know that Oskar gene exist")

(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
pathcsv = options.pathcsv
genomeSearch = options.genomeSearch


if genomeFolder == "None":
	print("Genome Folder must be provided")
	sys.exit(1)

if pathcsv == "None":
	print("curated_InsectGenomesInfos.csv file must be provided")
	sys.exit(1)

if genomeSearch == "None":
	print("GenomeSearch file must be provided")
	sys.exit(1)

pre_insectDataFrame, allGenomesInfos = getAllGenomesInfos(pathcsv, genomeFolder)
buildpipelineGuide(allGenomesInfos, pre_insectDataFrame)