#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-13-2017
Modified date:

Description: Add fasta and gff file of all refseq sequences
in folder for running Augustus training

'''

import os, sys, re
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-r", "--refseq_path", dest="refseqPath", default="None",
                  help="[Required] Location of the directory where are Refseq genomes")
parser.add_option("-o", "--output_path", dest="augustusFolderPath", default="None",
                  help="[Required] Location of the directory where will be built Augustus directories ")

(options, args) = parser.parse_args()

refseqPath = options.refseqPath
augustusFolderPath = options.augustusFolderPath

if refseqPath == "None":
	print("Refseq path must be provided.\n -h for more information")
	sys.exit(1)

if augustusFolderPath == "None":
	print("Refseq path must be provided.\n -h for more information")
	sys.exit(1)

genome_folder = os.listdir(refseqPath)

for folder in genome_folder :
    path = os.path.join(refseqPath, folder)

    file_list = os.listdir(path) :

    for genomeFile in file_list :
