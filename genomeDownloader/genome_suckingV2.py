#!/usr/bin/env python 

'''
Author: Savandara Besse and Leo Blondel

Created date: 03-14-2017
Modified date: 03-15-2017

Description: This script allows to download all genomes folders 
corresponding to accession numbers given in the file 
curated_insectGenomesInfos.csv

'''

import pandas as pd
import re
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--csv_file", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file containing all the genomes ID in the collumn named 'Accession number'")
parser.add_option("-o", "--genomePath", dest="genomePath", default="None",
                  help="[Required] Location of the directory where the genomes will be downloaded")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv
genomePath = options.genomePath

if pathcsv == "None":
	print("List of genome to dowload must be provided.\n -h for more information")
	sys.exit(1)

datatofetch = pd.read_csv(pathcsv)
gIDtofetch = datatofetch['Accession number'].values

if not os.path.isdir(os.path.join(genomePath, 'AllGenomes')):
	os.mkdir(os.path.join(genomePath, 'AllGenomes'))

for i in range(len(gIDtofetch)):
	arg = re.findall("(^[A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3})", gIDtofetch[i])[0]
	path = os.path.join(genomePath, 'AllGenomes', ''.join(arg))
	if not os.path.isdir(path):
		os.mkdir(path)
	os.system('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/ '.format(arg[0], arg[1], arg[2], arg[3])+path)