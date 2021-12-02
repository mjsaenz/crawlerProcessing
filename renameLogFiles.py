import sys
import  glob
import shutil
import os

path = "/home/c2i/gss_logs/"
outDirRoot = '/home/c2i/archive'
#####################################
flist = sorted(glob.glob(os.path.join(path, '*.gssbin')))

for fname in flist:
    date = os.path.basename(fname).split("_")[0]
    year = date[0:4]  # grab year from date
    outDir = os.path.join(outDirRoot, year, date)
    if not os.path.exists(outDir):
        os.makedirs(outDir)  # make the directory (including super directories)
    # now find all files that match the "date"
    moveFlist = glob.glob(os.path.join(os.path.dirname(fname), f"*{date}_*"))

    for f in moveFlist:
        try:
            shutil.move(f, outDir)
        except shutil.Error:
            #likely an error for "file already exists"
            print(f"{os.path.basename(fname)} exists in the target directory ")

