#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-13-2017
Modified date:

Description: Add fasta and gff file of all refseq sequences
in folder for running Augustus training

'''

import os, sys, re
import subprocess
from subprocess import Popen
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
    print('Building {} folder in progress ...'.format(folder))
    path = os.path.join(refseqPath, folder, folder+'_genomic.fna.gz')
    if os.path.isfile(path) :
        if not os.path.isfile(os.path.join(refseqPath, folder, folder+'_genomic.fna')):
            P = Popen(['gunzip',path])
            ret = P.wait()
            if ret != 0:
                print("Error Gunzipping !")

    path = os.path.join(refseqPath, folder, folder+'_genomic.fna')
    if os.path.isfile(path):
        if not os.path.isdir(os.path.join(augustusFolderPath, folder)):
            os.mkdir(os.path.join(augustusFolderPath, folder))
        os.system('cp {} {}'.format(path, os.path.join(augustusFolderPath, folder)))

    path = os.path.join(refseqPath, folder, folder+'_genomic.gff.gz')
    if os.path.isfile(path) :
        if not os.path.isfile(os.path.join(refseqPath, folder, folder+'_genomic.gff')):
            P = Popen(['gunzip',path])
            ret = P.wait()
            if ret != 0:
                print("Error Gunzipping !")

    path = os.path.join(refseqPath, folder, folder+'_genomic.gff')
    if os.path.isfile(path):
        path = os.path.join(refseqPath, folder, folder+'_genomic.gff')
        os.system('cp {} {}'.format(path, os.path.join(augustusFolderPath, folder)))
