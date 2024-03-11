import sys
import glob
import os
import os.path

year='2018'


files = glob.glob('*{0:s}'.format(year))
for file in files :
    #if 'ZH' not in file and 'HZ' not in file : continue
    print file
    #if '2018C' not in file : continue
    #if 'DY' not in file : continue
    if not os.path.isdir(file) : continue
    os.chdir('./{0:s}'.format(file))
    print("cwd={0:s}".format(os.getcwd())) 
    jdls = glob.glob('*.jdl')
    for jdl in jdls :
        ff=jdl.replace('jdl','root')
        #print("Command={0:s}".format(command))
        rootfilename = ff.replace('_{0:s}'.format(year),'')
        feos = rootfilename.replace('all_','')
        cf = os.path.isfile('/eos/uscms/store/user/lpcsusyhiggs/ntuples/Suu/{0:s}/{1:s}'.format(file,feos))
        print '/eos/uscms/store/user/lpcsusyhiggs/ntuples/Suu/{0:s}/{1:s}'.format(file,feos)

        #print '/eos/uscms/store/user/lpcsusyhiggs/ntuples/nAODv9/2Lep/{0:s}/{1:s}'.format(file,feos), 'feos--<', feos
        cfs = os.path.isfile('{0:s}.submitted'.format(jdl))
        #if  cf and cfs: print 'The .root exists:', ff, ' I wont submit'
        #else: 
        if not cf or not cfs:
            command = "condor_submit {0:s}".format(jdl)
            commandd = "touch {0:s}.submitted".format(jdl)
            print 'sending...', ff, command, jdl
            os.system(command)
            os.system(commandd)
    os.chdir('..')
    print("cwd={0:s}".format(os.getcwd()))

