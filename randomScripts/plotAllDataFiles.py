import pandas as pd
import glob, os, sys
from matplotlib import pyplot as plt

######################### plot all Files ##############################
outprefix = "singleFileAllCSVdump/plotsOut"
singlePlotThresh = 8
flist = glob.glob("singleFileAllCSVdump/*.csv")
plt.ioff()
for fname in flist:
    df = pd.read_csv(fname, header=4)
    base = fname.split('gssbin_')[-1].split('.')[0]

    badKeys= []
    for key in df.keys():  # plot
        print('file: {}\n   key: {}'.format(base, key))
        try:
            plt.figure()
            plt.title('{}: {}'.format(base, key))
            df[key].plot(marker='.')
            plt.tight_layout()
            plt.savefig(os.path.join(outprefix, base+"_"+key + '.png'))
            plt.close()
        except TypeError:
            badKeys.append(key)
            plt.close()
            continue
        finally:
            plt.figure()
            plt.title("FILE: {}".format(base), fontweight='bold', fontsize=12)
            for i in range(len(badKeys)):
                plt.text(0, i, "{}".format(badKeys[i]), fontweight='bold', fontsize=9)
                plt.text(1, i, "{}".format(df[badKeys[i]][:3]), fontsize=9)

            plt.ylim(0, len(badKeys))
            plt.xlim(0, 3)
            plt.axis('off')
            plt.tight_layout()
            plt.savefig(os.path.join(outprefix, base+"___BADKEYS.png"))
            plt.close()

            if df.keys().size - len(badKeys) < singlePlotThresh:
                df.plot(subplots=True, figsize=(10, 14), title=base)
                plt.savefig(os.path.join(outprefix, base + "___summary.png"))
                plt.close()

