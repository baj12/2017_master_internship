'''
Author: Savandara Besse

Created date: 03/28/2017
Modified date: 03/29/2017

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

    for index in df.index :
        if df['#EST'][index] != 0 :
            print('Downloading EST: '+ df['Species name'][index])
            search_handle = Entrez.esearch(db = "nucest", term = df['Species name'][index], field='Organism', retmax =df['#EST'][index])
            id_list = Entrez.read(search_handle)['IdList']
            search_handle.close()

            out_handle = open(estpath+'/'+df['Species name'][index]+'_est.fna',"w")
            for elem in range(len(id_list)) :
                print('Downloading '+str(elem+1) +'/'+str(len(id_list))+' in progress')
                net_handle = Entrez.efetch(db="nucest", id=id_list[elem], rettype="fasta", retmode="text")
                out_handle.write(net_handle.read())
            net_handle.close()
            out_handle.close()
            print(df['Species name'][index]+'_est.fna saved')

# Main

parser = OptionParser()
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file hmmOrder.csv")
parser.add_option("-e", "--estpath", dest="estpath", default="None",
                  help="[Required] Location of the directory where the est will be downloaded")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv
estpath = options.estpath

download_EST(pathcsv, estpath)
