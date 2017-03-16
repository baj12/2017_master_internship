#!/usr/bin/env python 

'''

Author: Savandara Besse and Leo Blondel

Created date: 03-14-2017
Modified date: 03-15-2017

Description: This script allows to download all genomes folders corresponding to accession numbers given in the file curated_insectGenomesInfos.csv
Notes: Don't forget to create your result 'AllGenomes' before running this script
'''

import pandas as pd
import re
import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--csv_file", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file containing all the genomes ID in the collumn named 'Accession number'")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv

if pathcsv == "None":
	print("List of genome to dowload must be provided.\n -h for more information")
	sys.exit(1)

datatofetch = pd.read_csv(pathcsv)
gIDtofetch = datatofetch['Accession number'].values

for i in range(len(gIDtofetch)):
	arg = re.findall("(^[A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3})", gIDtofetch[i])[0]
	os.system('mkdir ../AllGenomes/'+''.join(arg)+ '; rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}/{}/{}/'.format(arg[0], arg[1], arg[2], arg[3])+' ../AllGenomes/'+''.join(arg)+'/')