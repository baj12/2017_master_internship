#!/usr/bin/env python 

'''

Author: Savandara Besse

Created date: 03-16-2017
Modified date: 

Description: This script cleans the gff output file from Exonerate

'''

import sys
import os
import re
from optparse import OptionParser
import progressbar


def deleteUselessLines(allLines):
	for line in allLines[2:]:

		if re.match(r'^#.*',line) is not None :
			allLines.remove(line)

	return allLines

def createCuratedFile(finalLines):
	f = open("protein_alignments.gff3", 'w')
	f.write("##gff-version 3\n")

	for line in finalLines[2:] :
		f.write(line)
		f.write('\n')

	f.close()

# Main

parser = OptionParser()
parser.add_option("-p", "--gff_file", dest="pathgff", default="None",
                  help="[Required] Location of the gff file that we want to clean'")

(options, args) = parser.parse_args()

pathgff = options.pathgff

if pathgff == "None":
	print("Gff file must be provided.\n -h for more information")
	sys.exit(1)

f = open(pathgff)
allLines = [x.strip() for x in f.readlines()]
f.close()

finalLines = deleteUselessLines(allLines)
createCuratedFile(finalLines)
