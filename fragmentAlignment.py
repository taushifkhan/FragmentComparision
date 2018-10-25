
import sys, os
import subprocess
import re
import pickle
import numpy as np
# from tqdm import tqdm
tmDir = "/Users/taushif/Taushif_Ib2/EarlyFolding/pisces/src/alignmentTool"
cwd = os.getcwd()+'/'

def tmalign(fl1,fl2,absPath=0):
    """
    @param:
    fl1 - > file name/path of file1
    fl2 - > file name/path of file 2
    absPath -> if fl1 nad fl2 are absolute path . If not current work directory will be appended.

    @return:
    0: if files are not found
    scrs [dictionary] :  dictionary with rmsd, se_Id, TmScore
    """
    if absPath:
        f1path = fl1
        f2path = fl2
    else:
        f1path = cwd+fl1
        f2path = cwd+fl2

    if os.path.isfile(f1path) or os.path.isfile(f2path):
        cmd = "%s/TMalign %s %s"%(tmDir,f1path,f2path)
        proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        if err:
            return  0
        else:
            l = []
            for i in out.split("\n"):
                if re.search("=",i):
                    l.append(i)
            scrs = {}

            scrs[l[0].split(",")[0].split("=")[0].strip()] = int(l[0].split(",")[0].split("=")[1])
            scrs[l[0].split(",")[1].split("=")[0].strip()] = float(l[0].split(",")[1].split("=")[1])
            scrs[l[0].split(",")[2].split("=")[0].strip()] = float(l[0].split(",")[2].split("=")[2])
            scrs[l[1].split("=")[0].strip()] = float(l[1].split("=")[1].split()[0])
            return scrs
    else:
        return 0

def ramRMSD(frag1,frag2):
    """
    Distance between dihedral 
    @param
    frag1: List pf arrays for phi psi angles for fragment 1
    frag2: List of array for phi psi for fragment 2

    @return
    delta T
    Refernce :
    # 1.1. Karpen ME, 1989,proteins
    """
    d_= np.array(frag1)-np.array(frag2)
    d_mod = []
    for row in d_:
        tmp =[]
        for i in row:
            if i**2 > 180**2:
                tmp.append((360-abs(i))**2)
            else:
                tmp.append(i**2)
        d_mod.append(tmp)
        
    d_mod = np.array(d_mod)
    return np.sqrt(sum(np.sum(d_mod,axis=1))/float(2*len(frag1)))

def BCsearch_ASD(fl1,fl2):
    scrs = {'bc':999,'asd':0}
    cmd_Ca_fl1 = """grep ATOM %s| grep "CA " |cut -c30-55>%s.xyz"""%(fl1,fl1)
    proc1 = subprocess.Popen([cmd_Ca_fl1], stdout=subprocess.PIPE, shell=True)
    (out1,err1) = proc1.communicate()
    cmd_Ca_fl2 = """grep ATOM %s| grep "CA " |cut -c30-55>%s.xyz"""%(fl2,fl2)
    proc2 = subprocess.Popen([cmd_Ca_fl2], stdout=subprocess.PIPE, shell=True)
    (out2,err2) = proc2.communicate()
    bcSpath = "./alignmentTool/BCSearch"
    asdPath = "./alignmentTool/asd.py"
    if err1 or err2:
        print "error"

    else:
        bcS_cmd = "%s -if %s -ip %s"%(bcSpath,fl1+'.xyz',fl2+'.xyz')
        asd_cmd = "python %s %s %s"%(asdPath,fl1+'.xyz',fl2+'.xyz' )
        
        proc_bcS = subprocess.Popen([bcS_cmd], stdout=subprocess.PIPE, shell=True)
        proc_asd = subprocess.Popen([asd_cmd], stdout=subprocess.PIPE, shell=True)
        result = proc_bcS.communicate()
        result_asd = proc_asd.communicate()
        try:
            scrs['asd'] = float(result_asd[0].split(":")[1].strip())
        except:
            pass
        try:
            scrs['bc'] = float(result[0].split(" ")[2].strip())
        except:
            pass
        
        subprocess.Popen(["rm %s.xyz %s.xyz"%(fl1,fl2)], stdout=subprocess.PIPE, shell=True)
            
    return scrs


def computeAlignments(frag_names,logFl_name):
    fragmentAln = []
    err = []
    bcactual = []
    for i in tqdm(range(len(frag_names)-1)):
        fl1 = frag_names[i]
        fl1_n = fl1[11:].strip(".pdb")
        logFl = open(logFl_name,'a')
        logFl.write("Working on %s .. [to compute %d]\n"%(fl1, len(frag_names)-i))
        for j in range(i+1,len(frag_names)):
            fl2 = frag_names[j]
            if os.path.isfile(fl1) and os.path.isfile(fl1):
                fl2_n = fl2[11:].strip(".pdb")
                ascore = tmalign(fl1, fl2)
                bcS = BCsearch_ASD(fl1, fl2)

                fragmentAln.append([fl1_n,fl2_n,ascore[' RMSD'],ascore['TM-score']\
                                    ,ascore[' Seq_ID'],ascore['Aligned length'],bcS['bc'],bcS['asd']])
            else:
                next
            logFl.close()
    return fragmentAln

def get_fragments(fragDict):
    resultfiles= {}

    for k in fragDict.keys():
        if k >=8 and k <=12:
            logfl_name = 'log_%d.txt'%(k)

            aln_result = {}
            d = fragDict[k]
            dfile = "%d_alignresult_%d.pkl"%(k,len(d))
            aln_result = computeAlignments(d, logfl_name)
            pickle.dump(aln_result,open(dfile,"wb"))
            resultfiles[k]=dfile

    return resultfiles


def main():
    fdict = pickle.load(open("frag_pathDict.pkl","rb"))
    alnflPath = get_fragments(fdict)
    print "Files in :"
    for k in alnflPath.keys():
        print k,":", alnflPath[k]

if __name__ == "__main__":
    main()