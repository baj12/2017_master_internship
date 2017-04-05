'''
Author: Savandara Besse

Created date: 03/31/2017
Modified date:

Description: Build each folder for doing Maker analysis on the 13 HMM orders
For each folder, we have to get :
- the genome fasta file
- the TSA or mRNA fasta file
- the protein database fasta file
'''
import re
import os, sys
import pandas as pd
import subprocess
from subprocess import Popen
from optparse import OptionParser


def unzipFile(genomefolder, currentG, gz_suffix, suffix):
    genomeidFolder = os.listdir(genomefolder)

    toFind = re.findall(r'^(GCA)_([0-9]{9}).[0-9]{1}', currentG)[0]
    id_folder = toFind[0] + toFind[1]

    path = os.path.join(genomefolder, id_folder)
    if not os.path.isfile(path) :
        genomeVersion = os.listdir(path)

        for  genome in genomeVersion :
            if currentG in genome :
                print("Opening Genome: {}".format(genome)+gz_suffix)
                path = os.path.join(genomefolder, id_folder, genome, genome+gz_suffix)
                if os.path.isfile(path):
                    if not os.path.isfile(os.path.join(genomefolder, id_folder, genome, genome+suffix)):
                        P = Popen(['gunzip',path])
                        ret = P.wait()
                        if ret != 0:
                            print("Error Gunzipping !")
                path = os.path.join(genomefolder, id_folder, genome, genome+suffix)
                if os.path.isfile(path):
                    return path


def buildFolder(pathcsv, estpath, proteinpath):
    df = pd.read_csv(pathcsv)

    for index in df.index :
        nameFolder = df['Order name'][index]+'_maker'
        if not os.path.isdir(nameFolder):
            os.mkdir(nameFolder)

        estFolder = os.listdir(estpath)
        for est_file in estFolder :
            if df['Species name'][index] in est_file :
                est_file = est_file.replace(' ', '\ ')
                os.system('cp {}/{} ./{}/'.format(estpath,est_file,nameFolder))

        genomeID = df["Accession number"][index]
        genomepath = unzipFile(genomefolder, genomeID, "_genomic.fna.gz", "_genomic.fna")
        os.system('cp {} ./{}/'.format(genomepath,nameFolder))
        os.system('cp {} ./{}/'.format(proteinpath,nameFolder))

    print('Maker folder built')


# Main

parser = OptionParser()
parser.add_option("-p", "--pathcsv", dest="pathcsv", default="None",
                  help="[Required] Location of the csv file hmmOrder.csv")
parser.add_option("-g", "--genome_folder", dest="genomefolder", default="None",
                  help="[Required] Location of the folder containing all the genomes")
parser.add_option("-e", "--estpath", dest="estpath", default="None",
                  help="[Required] Location of the directory where are est")
parser.add_option("-d", "--proteinpath", dest="proteinpath", default="None",
                  help="[Required] Location of the directory where is protein database file")

(options, args) = parser.parse_args()

pathcsv = options.pathcsv
estpath = options.estpath
genomefolder = options.genomefolder
proteinpath = options.proteinpath

buildFolder(pathcsv, estpath, proteinpath)
