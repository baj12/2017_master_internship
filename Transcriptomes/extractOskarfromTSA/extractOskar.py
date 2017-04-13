#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-11-2017
Modified date: 04-12-2017

Description: Extract all protein sequences which are supposed
to have an oskar sequence

'''

import pandas as pd
import re, os, sys
import progressbar
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

parser = OptionParser()
parser.add_option("-i", "--protein_path", dest='protein_path', default="None",
                  help="[Required] Location of protein folder")
parser.add_option("-p", "--oskar_result_path", dest='oskar_result_path', default="None",
                  help="[Required] Location of oskar_result.csv")

(options, args) = parser.parse_args()
protein_path = options.protein_path
oskar_result_path = options.oskar_result_path

df = pd.read_csv(oskar_result_path)

proteinOskar = {}

for index in df.index :
    if df['oskar'][index] == 'yes':
        tsa_id = df['tsa'][index]
        proteinOskar[tsa_id] = {}
        tsa = re.findall(r'([A-Z]{4}[0-9]{8}.[0-9]{1}_[0-9]{1})', df['oskar_sequence'][index])
        seq_list = []
        for i in range(len(tsa)):
            seq_list.append(tsa[i])
        proteinOskar[tsa_id]['oskar'] = seq_list
        proteinOskar[tsa_id]['organism'] = df['organism'][index]

proteinPath = os.listdir(protein_path)

my_records = []
bar = progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ', ])
for protein_file in bar(proteinPath) :
    for key in sorted(proteinOskar.keys()):
        if key in protein_file :
            for sequence in proteinOskar[key]['oskar'] :
                for seq_record in SeqIO.parse('{}/{}'.format(protein_path,protein_file), 'fasta'):
                    tsa = re.findall(r'([A-Z]{4}[0-9]{8}.[0-9]{1}_[0-9]{1})', seq_record.id)[0]
                    if sequence == tsa :
                        tsa_id = re.findall(r'([A-Z]{4}[0-9]{8}.[0-9]{1})', sequence)[0]


                        tmp = SeqRecord( Seq(str(seq_record.seq), generic_protein), id='{}|transeq_tsa|'.format(tsa_id), description=str(proteinOskar[key]['organism']) )

                        my_records.append(tmp)

SeqIO.write(my_records, "tsa_oskar_sequence.faa", "fasta")
