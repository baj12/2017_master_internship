#!/usr/bin/env python

'''
Author: Savandara Besse and Leo Blondel

Created date: 04-05-2017
Modified date: 03-15-2017

Description: This script allows to download all genomes folders
corresponding to accession numbers given in the file
tsa_selector.csv

'''

import os, sys, re
import progressbar
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

if not os.path.isdir(os.path.join(genomePath, 'diptera_TSA')):
	os.mkdir(os.path.join(genomePath, 'diptera_TSA'))

f = open(pathcsv)
alldata = [x.strip() for x in f.readlines()]
f.close()

tsaDict = {}

for line in alldata[1:] :
    line = line.split('\t')
    tsaDict[line[2]] = line[3]


bar = progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ', ])
path = os.path.join(genomePath, 'diptera_TSA')
for key in bar(tsaDict.keys()) :
    init = re.findall("^([A-Z]{4})", key)[0]
    version = tsaDict[key]
    os.system('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genbank/tsa/tsa.{}.{}.fsa_nt.gz '.format(init, version)+path)
