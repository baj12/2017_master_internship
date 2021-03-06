#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-08-2017
Modified date : 04-12-2017

Description: This script compares result files from hmmsearch to identify
if an osk and an lotus is in the same sequence

'''

import sys, os, re
import json
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--oskar_search_path", dest='oskar_search_path', default="None",
                  help="[Required] Location of the hmmsearch result folder")

(options, args) = parser.parse_args()
oskar_search_path = options.oskar_search_path

oskar_search = os.listdir(oskar_search_path)

withOskar = 0
withoutOskar = 0

oskarDict = {}
for folder in oskar_search :
    path = os.path.join(oskar_search_path, folder)

    for file in os.listdir(path) :
        if "osk" in file :
            osk = file
        elif "lotus" in file :
            lotus = file

    osk_file = open('{}/{}/{}'.format(oskar_search_path, folder, osk))
    oskLines = [x.strip() for x in osk_file.readlines()]
    osk_file.close()

    lotus_file = open('{}/{}/{}'.format(oskar_search_path, folder, lotus))
    lotusLines = [x.strip() for x in lotus_file.readlines()]
    lotus_file.close()

    osk_idList = []
    for line in oskLines[3:] :
        if "#" not in line :
            orgn = re.findall(r'TSA: ([A-Za-z]* [a-z]*)',line)[0]
            line = line.split('-')
            osk_idList.append(re.findall(r'[A-Z]{4}[0-9]{8}.[0-9]{1}_[0-9]{1})',line[0]))

    lotus_idList = []
    for line in lotusLines[3:] :
        if "#" not in line :
            line = line.split('-')
            lotus_idList.append(re.findall(r'[A-Z]{4}[0-9]{8}.[0-9]{1}_[0-9]{1})',line[0]))

    tsa = re.findall(r'^([A-Z]{4}.[0-9]{1})',osk)[0]

    '''
    for line in alldata[1:] :
        line = line.split('\t')
        tsaDict[line[2]] = line[3]
    '''


    data = {}
    data['metadata'] = []
    data['metadata'].append( {
    })
    if len(osk_idList) == 0 and len(lotus_idList) == 0 :
        withoutOskar += 1
        oskarDict[tsa] = {}
        oskarDict[tsa]['orgn'] = orgn
        oskarDict[tsa]['id'] = ''
        oskarDict[tsa]['osk'] = 'no'
        oskarDict[tsa]['lotus'] = 'no'
        oskarDict[tsa]['oskar'] = 'no'


    elif len(osk_idList) != 0 and len(lotus_idList) == 0 :
        withoutOskar += 1
        oskarDict[tsa] = {}
        oskarDict[tsa]['orgn'] = orgn
        oskarDict[tsa]['id'] = ''
        oskarDict[tsa]['osk'] = 'yes'
        oskarDict[tsa]['lotus'] = 'no'
        oskarDict[tsa]['oskar'] = 'no'

    elif len(osk_idList) == 0 and len(lotus_idList) != 0 :
        withoutOskar += 1
        oskarDict[tsa] = {}
        oskarDict[tsa]['orgn'] = orgn
        oskarDict[tsa]['id'] = ''
        oskarDict[tsa]['osk'] = 'no'
        oskarDict[tsa]['lotus'] = 'yes'
        oskarDict[tsa]['oskar'] = 'no'

    else :
        withOskar += 1
        id_list = []
        for osk_id in osk_idList :
            if osk_id in lotus_idList :
                id_list.append(osk_id)
                oskarDict[tsa] = {}
                oskarDict[tsa]['orgn'] = orgn
                oskarDict[tsa]['id'] = id_list
                oskarDict[tsa]['osk'] = 'yes'
                oskarDict[tsa]['lotus'] = 'yes'
                oskarDict[tsa]['oskar'] = 'yes'
            else :
              oskarDict[tsa] = {}
              oskarDict[tsa]['orgn'] = orgn
              oskarDict[tsa]['id'] = ''
              oskarDict[tsa]['osk'] = 'yes'
              oskarDict[tsa]['lotus'] = 'yes'
              oskarDict[tsa]['oskar'] = 'no'

f = open('oskar_result.csv', 'w')
f.write('tsa,organism,osk,lotus,oskar,oskar_sequence\n')

for key in sorted(oskarDict.keys()):
    f.write('{},{},{},{},{},'.format(key, oskarDict[key]['orgn'], oskarDict[key]['osk'], oskarDict[key]['lotus'], oskarDict[key]['oskar']))
    for elem in oskarDict[key]['id'] :
        f.write(elem)
    f.write('\n')

f.close()
