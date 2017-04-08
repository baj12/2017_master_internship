#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-07-2017

Description: This script allows run hmmsearch

'''

import os, sys, re
import progressbar
from subprocess import Popen
import subprocess
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--proteinPath", dest="proteinPath", default="None",
                  help="[Required] Location of the protein folder")
parser.add_option("-o", "--oskModel", dest="oskModel", default="None",
                  help="[Required] Location of osk hmm")
parser.add_option("-l", "--lotusModel", dest="lotusModel", default="None",
                  help="[Required] Location of lotus hmm")

(options, args) = parser.parse_args()

proteinPath = options.proteinPath
oskModel = options.oskModel
lotusModel = options.lotusModel

protein_folder = os.listdir(proteinPath)

if not os.path.isdir('oskar_search'):
    os.mkdir('oskar_search')

bar = progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ', ])
for protein in bar(protein_folder) :
    if not os.path.isdir('./oskar_search/{}_oskar_search'.format(protein)):
    	os.mkdir('./oskar_search/{}_oskar_search'.format(protein))

    path = os.path.join(proteinPath, protein)
    if os.path.isfile(path) :
        output = re.findall(r'^tsa.([A-Z]{4}.[0-9]{1})',protein)[0]
        os.system('hmmsearch --tblout ./oskar_search/{}_oskar_search/{}_osk_search.txt {} {}'.format(protein,output,oskModel,protein))
        os.system('hmmsearch --tblout ./oskar_search/{}_oskar_search/{}_lotus_search.txt {} {}'.format(protein,output,lotusModel,protein))
