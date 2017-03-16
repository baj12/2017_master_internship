#!/usr/bin/env python

'''

Author : Savandara Besse

Created date : 03-08-2017
Modified date : 03-09/2017

Description : This script allows to create a csv file composed 
of all the insect genomes informations with these specific columns :
Species id - Species name, Families id, Families name, Order id, Order name, acession number

To update this table, you will need of the last version of this file eucaryotes: 
ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt

'''

import sys
import os
from optparse import OptionParser
from Bio import Entrez
import progressbar
Entrez.email = 'savandara.besse@gmail.com'

 	
parser = OptionParser()
parser.add_option("-p", "--csv_file", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file containing all the data of this link: 'ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt'")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv

if pathcsv == "None":
	print ("Maybe have you forget to download your data? \n -h for more information")
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
			insectDict[line[1]]["genome_id"] = line[8]

	return insectDict

def addData(insectDict) :
	bar = progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ', ])
	for key in bar(insectDict.keys()):
		tax_id = key

		handle = Entrez.efetch(db="taxonomy", id=tax_id)
		record = Entrez.read(handle) 

		for index in range(len(record)) : 
			lineage_ex = record[index]["LineageEx"] 
			for Dict in lineage_ex :

				if Dict["Rank"] == "family" :
					insectDict[tax_id]["family_name"] = Dict["ScientificName"]
					insectDict[tax_id]["family_id"] = Dict["TaxId"]

				if Dict["Rank"] == "order" : 
					insectDict[tax_id]["order_name"] = Dict["ScientificName"]
					insectDict[tax_id]["order_id"] = Dict["TaxId"]

			if "family_name" not in insectDict[tax_id]:
			    insectDict[tax_id]["family_name"] = "N/A"
			    insectDict[tax_id]["family_id"] = "N/A"

			if "order_name" not in insectDict[tax_id]:
			    insectDict[tax_id]["order_name"] = "N/A"
			    insectDict[tax_id]["order_id"] = "N/A"
					
		handle.close()

	return insectDict


def createTable(insectDict):

	f = open("insectGenomesInfos.csv", 'w')
	f.write("Species id,Species name,Family id,Family name,Order id,Order name,Accession number\n")

	id_list = ["sp_name","family_id", "family_name", "order_id", "order_name"]

	for key in insectDict.keys():
		f.write(str(key))
		f.write(",")

		for id in range(len(id_list)) :
			value = insectDict[key][id_list[id]]

			if type(value) is int : 
				f.write(str(value))
				f.write(",")

			else :
				f.write(value)
				f.write(",")

		f.write(insectDict[key]["genome_id"])
		f.write("\n")

	f.close()


# Main

insectDict = checkIfInsect(alldata)
addData(insectDict)
createTable(insectDict)
