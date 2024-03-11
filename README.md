Simple scripts to create Suu skims

## Introduction
** Step 1: Prepare a csv with the nAOD samples to be skimmed. Make sure that there are no empty lines, and that the following format is 

```TTTo2L2Nu, Top, 87.3, 1, 1,,/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM```


the first column is used for the creation of the folder on /eos area

** Step 2:execute 

> python makeMC.py -f file.csv -y 2018 -s Suu

that will create a .csh file that has to be executed and it will create folders with both the .sh and the jdl for condor submittion

don't forget to edit folders like 

```https://github.com/alkaloge/Suu/blob/master/Skimmer/makeCondor.py#L165C1-L165C79```

to make sure that files to be send with the job are under your local area

final .root file is written on eos and can be changed here

```https://github.com/alkaloge/Suu/blob/master/Skimmer/makeCondor.py#L144```

** Step 3: create jobs > . runMC_2018_Suu.csh

** Step 4: > cp subAllDirCondor.py Suu; cd Suu;

open the subAllDirCondor.py and control what should be sent. The script will make also a .submitted file for each job submitted on condor to avoid re-sending the job again. Should you need to resent the job, don't forget to delete the corresponding .submitted file
