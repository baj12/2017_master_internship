#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-21-2017
Modified date: 03-22-2017 
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

def buildDataFrame(pathcsv):
	dfi = pd.read_csv(pathcsv)
	df = dfi.sort_values(["Order id", "Family id"], ascending=[0,1])[['Order name', 'Family name', 'Species name', 'Accession number']]
	df = df.set_index('Species name')
	
	return df


def makeDecisions(infoDict, currentG):

	if infoDict[currentG]["Protein annotations"] == 'Yes' :
		infoDict[currentG]["TO DO"] = 'Find Oskar'
	elif infoDict[currentG]["GeneBank annotations"] == 'Yes':
		infoDict[currentG]["TO DO"] = 'Use Augustus'
	elif infoDict[currentG]["GeneBank annotations"] == 'No' and (infoDict[currenG]["GFF3 annotations"] == 'Yes' and infoDict[currenG]["Genome Fasta file"] == 'Yes'):
		infoDict[currentG]["TO DO"] = 'Convert into GeneBank file and use Augustus!'
	elif infoDict[currentG]["GFF3 annotations"] == 'No' and infoDict[currenG]["Genomic fasta file"] == 'Yes':
		infoDict[currentG]["TO DO"] = 'Use GeneMark + Exonerate + Evidence Modeller + Augustus'
	elif infoDict[currentG]["Genomic fasta file"] == 'No' :
		infoDict[currentG]["TO DO"] = 'Should look à the folder '+ currentG
	else : 
		infoDict[currentG]["TO DO"] = "Good question ..."

	return infoDict


def getGenomeInfos(genomeFolder, folder, currentG, GCA_number, infoDict):
	infoDict[currentG] = {}
	infoDict[currentG]['Accession number'] = GCA_number
			
	path = os.path.join(genomeFolder, folder, currentG)
	allInfos = os.listdir(path)

	for info in allInfos :
		if '_protein.faa.gz' in info :
			infoDict[currentG]['Protein annotations'] = 'Yes'
		if '_genomic.gbff.gz' in info :
			infoDict[currentG]['GeneBank annotations'] = 'Yes'
		if '_genomic.gff.gz' in info :
			infoDict[currentG]['GFF3 annotations'] = 'Yes'
		if '_genomic.fna.gz' in info :
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

	


def getAllGenomesInfos(genomeFolder, df) :
	allGenomesInfos = []

	genomeidFolder = os.listdir(genomeFolder)
	for folder in [genomeidFolder[0],genomeidFolder[42],genomeidFolder[200]] :
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

			for genome in df['Accession number'] :
				regExpr = re.findall(r'^([A-Z]{3}_[0-9]{9}.[0-9]{1})', genome)[0]
				if hasToBeFound == regExpr :
					allGenomesInfos.append(getGenomeInfos(genomeFolder, folder, current, hasToBeFound, {}))

	return allGenomesInfos


def buildnewDataframe(infoDict, df)

# Main

parser = OptionParser()
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file curated_InsectGenomesInfos.csv")
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv
genomeFolder = options.genomeFolder

insectDataFrame = buildDataFrame(pathcsv)
result = getAllGenomesInfos(genomeFolder, insectDataFrame)
print(result)
