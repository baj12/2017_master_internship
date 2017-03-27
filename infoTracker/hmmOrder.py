#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-26-2017
Modified date: 03-27-2017

Description:
Create a csv file where can see for each organism and tell us :
- its order and its family
- if they have protein annotations file, its protein number
- if they have est in NCBI, its EST number

'''

import os, sys
import pandas as pd
import re
from optparse import OptionParser
from subprocess import Popen
import subprocess
from Bio import Entrez
Entrez.email = 'savandara.besse@gmail.com'


def unzipProteinFile(genomeFolder, currentG):
    genomeidFolder = os.listdir(genomeFolder)

    toFind = re.findall(r'^(GCA)_([0-9]{9}).[0-9]{1}', currentG)[0]
    id_folder = toFind[0] + toFind[1]

    path = os.path.join(genomeFolder, id_folder)
    if not os.path.isfile(path) :
        genomeVersion = os.listdir(path)

        for  genome in genomeVersion :
            if currentG in genome :
                print("Opening Genome: {}".format(genome))
                path = os.path.join(genomeFolder, id_folder, genome, genome+"_protein.faa.gz")
                if os.path.isfile(path):
                    if not os.path.isfile(os.path.join(genomeFolder, id_folder, genome, genome+"_protein.faa")):
                        P = Popen(['gunzip',path])
                        ret = P.wait()
                        if ret != 0:
                            print("Error Gunzipping !")
                path = os.path.join(genomeFolder, id_folder, genome, genome+"_protein.faa")
                if os.path.isfile(path):
                    return path

def countProtein(pathcsv):
    protInfo_list = []
    df = pd.read_csv(pathcsv)

    for species in df.index :
        genomeID = df["Accession number"][species]
        protDict = {}
        protDict[genomeID] = {}
        protDict[genomeID]["Accession number"] = genomeID

        if df["Protein annotations"][species] == "Yes" :
            path = unzipProteinFile(genomeFolder, genomeID)
            protein_nb = os.popen("grep '>' "+path+ " | wc -l").read().split('\n')[0]
            protDict[genomeID]["#Protein"] = protein_nb

        if df["Protein annotations"][species] == "No" :
            protDict[genomeID]["#Protein"] = "N/A"

        protInfo_list.append(protDict)

    prot_df = []

    for elem in protInfo_list :
        for key in elem.keys():
            prot_df.append(elem[key])

    df = df[['Order name', 'Family name', 'Species name', 'Accession number']]
    df_m1 = pd.DataFrame(prot_df)
    df_m2 = pd.merge(df, df_m1, how='inner', on=['Accession number'], sort=False, suffixes=('_x', '_y'), copy=True, indicator=False)

    return df_m2


def countEST(df):
    ESTInfo_list = []

    for species in df.index :
        orgn = df["Species name"][species]
        print("Looking for "+orgn+ " EST")
        handle = Entrez.esearch(db="nucest", term=orgn)
        record = Entrez.read(handle)["Count"]
        handle.close()

        estDict = {}
        estDict[orgn] = {}
        estDict[orgn]["Species name"] = orgn
        estDict[orgn]["#EST"] = record

        ESTInfo_list.append(estDict)

    est_df = []

    for elem in ESTInfo_list :
        for key in elem.keys():
            est_df.append(elem[key])

    df_m3 = pd.DataFrame(est_df)
    df_f = pd.merge(df, df_m3, how='inner', on=['Species name'], sort=False, suffixes=('_x', '_y'), copy=True, indicator=False)

    df_f = df_f[["Species name","Order name","Family name","#Protein","#EST"]]
    df_f = df_f.set_index('Species name')
    df_f.to_csv("prot_est_Infos.csv", sep=',', encoding="utf-8")

    return df_f

# Main

parser = OptionParser()
parser.add_option("-g", "--genome_folder", dest="genomeFolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file curated_InsectGenomesInfos.csv")

(options, args) = parser.parse_args()

genomeFolder = options.genomeFolder
pathcsv = options.pathcsv


if genomeFolder == "None":
	print("Genome Folder must be provided")
	sys.exit(1)

if pathcsv == "None":
	print("pipelineGuide file must be provided")
	sys.exit(1)

df = countProtein(pathcsv)
dataframe = countEST(df)
