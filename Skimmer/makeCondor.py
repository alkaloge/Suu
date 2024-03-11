from __future__ import division
import os
import subprocess
from decimal import *
import sys
getcontext().prec = 2

def getArgs() :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v","--verbose",default=0,type=int,help="Print level.")
    defDS = '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM '
    parser.add_argument("--dataSet",default=defDS,help="Data set name.") 
    parser.add_argument("--nickName",default='MCpileup',help="Data set nick name.") 
    parser.add_argument("-m","--mode",default='anaXRD',help="Mode (script to run).")
    parser.add_argument("-y","--year",default=2017,type=str,help="Data taking period, 2016, 2017 or 2018")
    parser.add_argument("-c","--concatenate",default=5,type=int,help="On how many files to run on each job")
    parser.add_argument("-s","--selection",default='ZH',type=str,help="select ZH or AZH")
    parser.add_argument("-j","--doSystematics",default='yes',type=str,help="do JME systematics")
    parser.add_argument("-l","--islocal",default='no',type=str,help="get list from /eos/ not DAS")
    parser.add_argument("-w","--weightsOnly",default='no',type=str,help="get list from /eos/ not DAS")
    parser.add_argument("-d","--isData",default=0,type=bool,help="data or mc")
    return parser.parse_args()

def beginBatchScript(baseFileName, Systematics) :
    outLines = ['#!/bin/sh\n']
    outLines.append("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
    outLines.append("setenv SCRAM_ARCH slc7_amd64_gcc820\n")
    outLines.append("eval `scramv1 project CMSSW CMSSW_10_6_5`\n")
    outLines.append("cd CMSSW_10_6_5/src\n")
    outLines.append("eval `scramv1 runtime -sh`\n")
    #outLines.append("tar -zxvf ${_CONDOR_SCRATCH_DIR}/taupog.tar.gz\n")
    outLines.append("scram b -j 4\n")
    #if Systematics :
    #if str(args.islocal.lower()=='yes') : 
    #    outLines.append("git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools\n")
    #    outLines.append("cp ${_CONDOR_SCRATCH_DIR}/branchselection.py PhysicsTools/NanoAODTools/python/postprocessing/framework/.\n")
    #    outLines.append("cp ${_CONDOR_SCRATCH_DIR}/keep_and_drop.txt PhysicsTools/NanoAODTools/python/postprocessing/framework/.\n") 
    #    outLines.append("cd PhysicsTools/NanoAODTools\n")
    #    outLines.append("scram b -j 4\n")
    outLines.append("echo ${_CONDOR_SCRATCH_DIR}\n")
    outLines.append("cd ${_CONDOR_SCRATCH_DIR}/CMSSW_10_6_5/src/\n")
    outLines.append("cp ${_CONDOR_SCRATCH_DIR}/* .\n")
    outLines.append("ls -altrh\n")
    outLines.append("echo 'this is the working dir' ${_CONDOR_SCRATCH_DIR}\n")
    return outLines

def getFileName(line) :
    tmp = line.split()[0].strip(',')
    fileName = tmp.strip()
    return fileName


args = getArgs()
era = str(args.year)
isdata = int(args.isData)

isdata = False
if 'Run' in str(args.nickName) : isdata = True

doJME  = args.doSystematics.lower() == 'true' or args.doSystematics.lower() == 'yes' or args.doSystematics == '1'

doWeightsOnly = args.weightsOnly.lower() =='yes' or args.weightsOnly.lower()=='1' or  args.weightsOnly.lower() == 'true'

period="B"
if 'Run2016' in args.dataSet or 'Run2017' in args.dataSet or 'Run2018' in args.dataSet: 
    poss = args.dataSet.find("Run")
    period = args.dataSet[int(poss)+7:int(poss)+8]
    print 'will set up', poss, period

#if 'HIPM' in args.dataSet: period='_preVFP'+period
if 'HIPM' in args.dataSet or '_pre' in args.dataSet: era=era+'_preVFP'


# sample query 
# dasgoclient --query="file dataset=/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8*/*/NANOAOD*" --limit=0   


query = '"file dataset={0:s}"'.format(args.dataSet)
if "USER" in str(args.dataSet) : query = '"file dataset={0:s}"'.format(args.dataSet+" instance=prod/phys03")
command = "dasgoclient --query={0:s} --limit=0  > fileList.txt".format(query)


if str(args.islocal.lower())=='yes' : command ='ls  /eos/uscms/store/group/lpcsusyhiggs/ntuples/AZh/nAODv9/ZH_JECs_{0:s}/CRAB_PrivateMC/{1:s}_{0:s}/*/*/*root   > fileList.txt'.format(str(args.year), args.nickName) 
#if str(args.islocal.lower())=='yes' or str(args.islocal.lower())=='1': command ='ls  /eos/uscms/store/group/lpcsusyhiggs/ntuples/nAODv7/JEC_{0:s}/*/{1:s}_{0:s}/*/*/*root   > fileList.txt'.format(str(args.year), args.nickName)

print("Running in {0:s} mode.  Command={1:s}".format(args.mode,command))
os.system(command)
    
files = open('fileList.txt','r').readlines()
if len(files) < 1 :
    print("***In makeCondor.py: Empty fileList.txt")
    exit()

scriptList = [] 
file=[]
dataset=[]

mjobs=args.concatenate


for nFiles, file in enumerate(files) :
     
    fileName=getFileName(file)
    if str(args.islocal.lower())=='yes' : fileName = fileName.replace('/eos/uscms','')
    if '#' not in fileName :  dataset.append(fileName)


counter=0

#executable = str(args.selection)
executable='SystWeights'


#doWeights=2
#mparts=2
mparts=mjobs

#if not doJME : mparts=1
wdir='/uscms_data/d3/alkaloge/MetStudies/nAOD/CMSSW_10_6_5/src/MakeSyst/test/'
ff = open("{0:s}/template_".format(wdir), "r")
textf = ff.read()
infile = "inFile.root"

#if str(args.islocal.lower())=='yes': 
if True:

    for nFile in range(0, len(dataset)) :
	#print("nFile={0:d} file[:80]={1:s}".format(nFile,file[:80]))
	fileName = getFileName(file)
	fileloop=dataset[nFile]
	outFileName = "{0:s}_{1:03d}.root".format(args.nickName,nFile+1)
        scriptName = "skim_{0:s}_file{1:s}.sh".format(args.nickName, str(nFile+1))
        runLines=[]
        runLines.append(textf)
       
	runLines.append("xrdcp root://cmsxrootd.fnal.gov/{0:s} inFile.root\n".format(fileloop)) 
	runLines.append("if [ ! -f inFile.root ] ; \n then \n xrdcp root://cms-xrd-global.cern.ch/{0:s} inFile.root\n fi \n".format(fileloop)) 
	extra=''
	if 'Run' in fileloop : runLines.append("python skim_template.py  0\n")
	else: runLines.append("python skim_template.py  1\n")
	runLines.append("rm inFile.root \n")

	runLines.append("xrdcp selected_events.root  root://cmseos.fnal.gov//store/group/lpcsusyhiggs/ntuples/Suu//{0:s}_{3:s}/skim_{0:s}_file{1:s}.root\n".format(args.nickName, str(nFile+1), outFileName, args.year)) 
	runLines.append("rm selected_events.root \n")
	runLines.append(" \n")

	runLines.append("rm inFile*.root\n")

	runLines.append("rm *.pyc\nrm *.so\nrm *.pcm\nrm *cc.d\n")
	runLines.append("rm *.ntup *.weights *.so\nrm *.pcm\nrm *cc.d\n")


	open(scriptName,'w').writelines(runLines)

	scriptList.append(scriptName)
	counter += mjobs



# now that .csh files have been generated make a list of corresponding .jdl files

#dir = '/uscms_data/d3/alkaloge/ZH/CMSSW_10_2_9/src/MC/'

wdir='/uscms_data/d3/alkaloge/MetStudies/nAOD/CMSSW_10_6_5/src/MakeSyst/test/'

dirMC = wdir+"/MC/"
dirZH = wdir+"/ZH/"
dirData = wdir+"/data/"
funcsDir = wdir+"/funcs/"
SVFitDir = wdir+"/SVFitDir/"
toolsDir = wdir+"/tools/"
pileupDir = wdir+"/pileup/"


print("dir={0:s}".format(dir))

for file in scriptList :
    base = file[:-3] 
    runLines = ['universe = vanilla\n']
    runLines.append('Executable = {0:s}\n'.format(file))
    runLines.append('Output = {0:s}.out\n'.format(base))
    runLines.append('Error = {0:s}.err\n'.format(base))
    runLines.append('Log = {0:s}.log\n'.format(base))
    #runLines.append('transfer_input_files = {0:s}ZH.py, {0:s}MC_{1:s}.root, {0:s}data_pileup_{1:s}.root, {0:s}MCsamples_{1:s}.csv, {0:s}ScaleFactor.py, {0:s}SFs.tar.gz, {0:s}cuts_{2:s}.yaml, '.format(dir,args.year, args.selection))
    #runLines.append('transfer_input_files = {0:s}{1:s}.py, '.format(dirZH,executable))
    #runLines.append('{0:s}MC_{1:s}.root, {0:s}data_pileup_{1:s}.root, '.format(dirMC,args.year, args.selection, executable))
    #runLines.append('{0:s}*txt, '.format(dirData))
    #runLines.append('{0:s}make_jmeMuons.py, {0:s}make_jmeElectrons.py, {0:s}make_jmeGjets.py, {0:s}keep_and_dropGjets.txt, {0:s}keep_and_drop.txt\n'.format(toolsDir))
    runLines.append('transfer_input_files ={0:s}/skim_template.py\n'.format(wdir))
    runLines.append('priority = 2\n')
    runLines.append('should_transfer_files = YES\n')
    runLines.append('when_to_transfer_output = ON_EXIT\n')
    runLines.append('x509userproxy = $ENV(X509_USER_PROXY)\n')
    runLines.append('Queue 1\n')
    open('{0:s}.jdl'.format(base),'w').writelines(runLines)


    
    
