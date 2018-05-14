
import os
import sys

batchDirName = "batch"

if len(sys.argv)>1:
    dirNames = sys.argv[1:]
else:
    dirNames = [batchDirName]

for dirname in dirNames:
    if dirname==batchDirName:
        dirpath = batchDirName
    else:
        dirpath = batchDirName+"/"+dirname
    listfiles = sorted(os.listdir(dirpath))
    for f in listfiles:
        if f.endswith(".json"):
            fpath = dirpath+"/"+f
            try:
                strcmd = """python3 run_optimization.py {}""".format(fpath)
                print("run batch:",strcmd)
                print(os.popen(strcmd).read())
            except:
                print("error in run, continue.",sys.exc_info()[0])

