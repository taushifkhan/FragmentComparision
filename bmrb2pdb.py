import urllib, json
from tqdm import tqdm
import xml.etree.ElementTree as ET

"""
link to API documentation: https://github.com/uwbmrb/BMRB-API

matchType:
BMRB Entry Tracking System: The entry is an exact match as tracked by the BMRB entry tracking system. There is a one-to-one correspondence between this queried entry and the provided BMRB ID.
BLAST Match - The entry was found during a routine BLAST search. It is similar to the queried entry in sequence but no other correlation is implied.
"""

def getPDB(bmrbId,matchType="BMRB Entry Tracking System"):
    linkAdd = 'http://webapi.bmrb.wisc.edu/v2/search/get_pdb_ids_from_bmrb_id/%s'%(bmrbId)
    f = urllib.urlopen(linkAdd)
    fJ = json.loads(f.read())

    for i in fJ:
        if i['match_type'] == matchType:
            return i['pdb_id']
        else:
            next
    # print "No one to one PBD id found for bmrb : ",bmrbId
    return 0

def getFasta(pdbId,chain='A'):
    url = 'https://www.rcsb.org/pdb/rest/das/pdbchainfeatures/sequence?segment=%s.%s'%(pdbId,chain)
    f = urllib.urlopen(url)
    fJ = ET.parse(f)
    r = fJ.getroot()
    for child in r:
        if child.attrib['moltype'] == 'Protein':
            fastaSeq = ">%s\n%s"%(child.attrib['id'],child.text)
            return fastaSeq
        else:
            next
    return 0



def main():
    import numpy
    dataDir = './'#'/Users/taushif/Taushif_Ib2/gmXsimulation/selectProteins/'
    bmrbFl = open(dataDir+'bmrbGoodEntries.txt','r').readlines()
    pdbFile = dataDir+'bmrbGoodEntries_pdb.txt'
    errFile = dataDir+'bmrbGoodEntries_error.txt'
    
    bmrb_pdb = []
    error = []
    flSize = len(bmrbFl)
    for k in tqdm(range(flSize)):
        i = bmrbFl[k]
        bid =  i.strip().split("_")[0][3:]
        pdb =  getPDB(bid)
        if pdb:
            bmrb_pdb.append([bid,pdb])
        else:
            error.append(bid)

    print "correct match:", len(bmrb_pdb), "No pdb:", len(error)
    numpy.savetxt(pdbFile,bmrb_pdb, fmt='%s')
    numpy.savetxt(errFile,error,fmt='%s')
    print "result in %s\nerror in :%s\n"%(pdbFile,errFile)

if __name__ == '__main__':
    import sys
    getFasta(sys.argv[1])
    # bmrbId = sys.argv[1] # example 7382
    # print bmrbId, getPDB(bmrbId)
    # main()


