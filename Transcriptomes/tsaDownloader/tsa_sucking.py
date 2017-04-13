#!/usr/bin/env python

'''
Author: Savandara Besse and Leo Blondel

Created date: 04-05-2017
Modified date: 04-12-2017

Description: This script allows to download all tsa fasta files
corresponding to accession numbers given in the file
tsa_selector.csv

'''

import os, sys, re
import progressbar
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--csvFile", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file containing all the genomes ID in the collumn named 'Accession number'")
parser.add_option("-o", "--tsaPath", dest="tsaPath", default="None",
                  help="[Required] Location of the directory where the genomes will be downloaded")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv
tsaPath = options.tsaPath

if pathcsv == "None":
	print("csv file must be provided.\n -h for more information")
	sys.exit(1)

if tsaPath == "None":
	print("tsa path must be provided.\n -h for more information")
	sys.exit(1)

if not os.path.isdir(os.path.join(tsaPath, 'tsa')):
	os.mkdir(os.path.join(tsaPath, 'tsa'))

f = open(pathcsv)
alldata = [x.strip() for x in f.readlines()]
f.close()

tsaDict = {}

for line in alldata[1:] :
    line = line.split('\t')
    tsaDict[line[2]] = line[3]


bar = progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ', ])
path = os.path.join(tsaPath, 'tsa')
for key in bar(tsaDict.keys()) :
    init = re.findall("^([A-Z]{4})", key)[0]
    os.system('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genbank/tsa/tsa.{}.1.fsa_nt.gz '.format(init)+path)
