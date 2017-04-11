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
	print("List of genome to download must be provided.\n -h for more information")
	sys.exit(1)

dataToFetch = pd.read_csv('GCA_GCF_insectGenomes.csv')

organismToFetch = []
for index in dataToFetch.index :
    if 'GCF' in dataToFetch['genome_id'][index] :
        organismToFetch.append(dataToFetch['sp_name'][index])

if not os.path.isdir(os.path.join(genomePath, 'Refseq_Genomes')):
	os.mkdir(os.path.join(genomePath, 'Refseq_Genomes'))

path = (os.path.join(genomePath, 'Refseq_Genomes'))
for organism in organismToFetch :
    organism = organism.split(' ')
    os.system('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/{}_{}/latest_assembly_versions/ {}'.format(organism[0], organism[1], path))
