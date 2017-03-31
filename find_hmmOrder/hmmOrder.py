#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-26-2017
Modified date: 03-30-2017

Description:
Create a csv file where can see for each organism and tell us :
- its order and its family
- if they have protein annotations file, its protein number
- if they have est in NCBI, its EST number
- if they have GFF3 annotations, its gene number

'''

import os, sys
import pandas as pd
import re
from optparse import OptionParser
from subprocess import Popen
import subprocess
from Bio import Entrez
Entrez.email = 'savandara.besse@gmail.com'


def unzipFile(genomeFolder, currentG, gz_suffix, suffix):
    genomeidFolder = os.listdir(genomeFolder)

    toFind = re.findall(r'^(GCA)_([0-9]{9}).[0-9]{1}', currentG)[0]
    id_folder = toFind[0] + toFind[1]

    path = os.path.join(genomeFolder, id_folder)
    if not os.path.isfile(path) :
        genomeVersion = os.listdir(path)

        for  genome in genomeVersion :
            if currentG in genome :
                print("Opening Genome: {}".format(genome)+gz_suffix)
                path = os.path.join(genomeFolder, id_folder, genome, genome+gz_suffix)
                if os.path.isfile(path):
                    if not os.path.isfile(os.path.join(genomeFolder, id_folder, genome, genome+suffix)):
                        P = Popen(['gunzip',path])
                        ret = P.wait()
                        if ret != 0:
                            print("Error Gunzipping !")
                path = os.path.join(genomeFolder, id_folder, genome, genome+suffix)
                if os.path.isfile(path):
                    return path


def countProtein(pathcsv):
    protInfo_list = []
    df_i = pd.read_csv(pathcsv)

    for species in df_i.index :
        genomeID = df_i["Accession number"][species]
        protDict = {}
        protDict[genomeID] = {}
        protDict[genomeID]["Accession number"] = genomeID

        if df_i["Protein annotations"][species] == "Yes" :
            path = unzipFile(genomeFolder, genomeID,"_protein.faa.gz","_protein.faa")
            protein_nb = os.popen("grep '>' "+path+ " | wc -l").read().split('\n')[0]
            protDict[genomeID]["#Protein"] = protein_nb

        if df_i["Protein annotations"][species] == "No" :
            protDict[genomeID]["#Protein"] = "N/A"

        protInfo_list.append(protDict)

    prot_df = []

    for elem in protInfo_list :
        for key in elem.keys():
            prot_df.append(elem[key])

    df_i = df_i[['Order name', 'Family name', 'Species name', 'Accession number']]
    df_m1 = pd.DataFrame(prot_df)
    prot_df = pd.merge(df_i, df_m1, how='inner', on=['Accession number'], sort=False, suffixes=('_x', '_y'), copy=True, indicator=False)

    return prot_df


def countGene(df, pathcsv):
    geneInfo_list = []
    df_i = pd.read_csv(pathcsv)

    for species in df_i.index :
        genomeID = df_i["Accession number"][species]
        protDict = {}
        protDict[genomeID] = {}
        protDict[genomeID]["Accession number"] = genomeID

        if df_i["GFF3 annotations"][species] == "Yes" :
            path = unzipFile(genomeFolder, genomeID, "_genomic.gff.gz", "_genomic.gff")
            protein_nb = os.popen("grep '##sequence-region' "+path+ " | wc -l").read().split('\n')[0]
            protDict[genomeID]["#Gene"] = protein_nb

        if df_i["GFF3 annotations"][species] == "No" :
            protDict[genomeID]["#Gene"] = "N/A"

        geneInfo_list.append(protDict)

    gene_df = []

    for elem in geneInfo_list :
        for key in elem.keys():
            gene_df.append(elem[key])

    df_m1 = pd.DataFrame(gene_df)
    prot_gene_df = pd.merge(df, df_m1, how='inner', on=['Accession number'], sort=False, suffixes=('_x', '_y'), copy=True, indicator=False)

    return prot_gene_df


def countEST(df):
    ESTInfo_list = []

    for species in df.index :
        orgn = df["Species name"][species]
        print("Looking for "+orgn+ " EST")
        handle = Entrez.esearch(db="nucest", term=orgn, field='Organism')
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

    df_f = df_f[["Species name","Order name","Family name",'Accession number',"#Gene","#Protein","#EST"]]
    df_f = df_f.set_index('Species name')
    df_f.to_csv("gene_prot_est_Infos.csv", sep=',', encoding="utf-8")

    print('gene_prot_est_Infos.csv done')


def build_hmm():
    df = pd.read_csv('gene_prot_est_Infos.csv')
    g = df.groupby(['Order name'])

    maxEST = g['#EST'].max().to_dict()

    f = open('hmmOrder.csv','w')
    f.write('Order name,Species name,Accession number,#EST\n')
    for index in df.index :
        for key in maxEST.keys():
            if df['Order name'][index] == key and df['#EST'][index] == maxEST[key]:
                f.write(df['Order name'][index])
                f.write(',')
                f.write(df['Species name'][index])
                f.write(',')
                f.write(df['Accession number'][index])
                f.write(',')
                f.write(str(df['#EST'][index])+'\n')
    f.close()

    print('hmmOrder.csv done')


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

prot_df = countProtein(pathcsv)
prot_gene_df = countGene(prot_df, pathcsv)
countEST(prot_gene_df)
build_hmm()
