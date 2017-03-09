#!/usr/bin/env python

'''

Author : Savandara Besse

Created date : 03-08-2017
Modified date : None

Description : This script allows to create a csv file composed 
of all the insect genomes informations with these specific columns :
Species names - Species id, Families names, Families id, Order names, Order id, genome id

To update this table, you will need of the last version of this file eucaryotes: 
ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt

'''

import pandas as pd
import os
from optparse import OptionParser
from Bio import Entrez

parser = OptionParser()
parser.add_option("-p", "--csv_file", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file containing all the data of this link: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt'")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv

if pathcsv == "None":
	print "Maybe have you forget to download your data? \n -h for more information"
	sys.exit(1)

f = open(pathcsv)
alldata = [x.strip() for x in f.readlines()]
f.close()

def checkIfInsect(alldata):
	insectDict = {}

	for line in alldata :
		line = line.split('\t')

		if line[5] == "Insects" :
			insectDict[line[1]] = {}
			insectDict[line[1]]["sp_name"] = line[0]
			insectDict[line[1]]["number_acession"] = line[8]

	return insectDict

def addData(insectDict) :
	for key in insectDict.keys():
		tax_id = key

		handle = Entrez.esearch(db="taxnomy", term=tax_id)
		record = Entrez.read(handle)

		record["Lineage"]


'''
Main
'''
checkIfInsect(alldata)

'''

def createTable(allDisease, diseaseDict):

	interestList = []

	for line in allDisease :
		line = line.split('\t')
		
		if line[0] == "Number Sign" :
			interestList.append(line)
		elif line[0] == "NULL" :
			interestList.append(line)
		elif line[0] == "Plus" :
			interestList.append(line)
		elif line[0] == "Percent" :
			interestList.append(line)
		elif line[0] == "Caret" :
			interestList.append(line)


	for disease in interestList :
		if diseaseDict.has_key(disease[1]): 
			diseaseDict[disease[1]]["name"] = disease[2]
			if disease[2] == "REMOVED FROM DATABASE" :
				del diseaseDict[disease[1]]
			move = re.search(r"MOVED TO (?P<id_number>\d+)", str(disease[2]))
			if move is not None : 
				for newDisease in interestList :
					if newDisease[1] == move.group('id_number') :
						diseaseDict[disease[1]]["name"] = newDisease[2]


	f = open("maladies.csv", 'w')
	f.write("OMIM number\t Disease Name\tTissue\n")

	for key in sorted(diseaseDict.keys()):
		f.write(str(key))
		f.write("\t")
		f.write(diseaseDict[key]["name"])
		f.write("\t")
		f.write(diseaseDict[key]["tissue"])
		f.write("\n")

	f.close()

'''