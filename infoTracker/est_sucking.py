'''
Author: Savandara Besse

Created date: 03-28-2017
Modified date:

Description: This script allows to download all est fasta for each organism
which choose as HMM for Maker pipeline
'''

import os
import pandas as pd
from optparse import OptionParser
from Bio import Entrez
Entrez.email = 'savandara.besse@gmail.com'

def download_EST(pathcsv, estpath):
    df = pd.read_csv(pathcsv)
    if df['#EST'][0] != 0 :
        print('Downloading EST: ')
        search_handle = Entrez.esearch(db = "nucest", term = df['Species name'][0],retmax =df['#EST'][0])
        id_list = Entrez.read(search_handle)['IdList']
        print(id_list)
        search_handle.close()

        out_handle = open(estpath+'/'+df['Species name'][0]+'_est.fna',"w")
        for elem in id_list :
            net_handle = Entrez.efetch(db="nucest", id=elem, rettype="fasta", retmode="text")
            out_handle.write(net_handle.read())
            net_handle.close()
        out_handle.close()
        print(df['Species name'][0]+'_est.fna saved')


# Main

parser = OptionParser()
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file curated_InsectGenomesInfos.csv")
parser.add_option("-e", "--estpath", dest="estpath", default="None",
                  help="[Required] Location of the directory where the est will be downloaded")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv
estpath = options.estpath

download_EST(pathcsv, estpath)
