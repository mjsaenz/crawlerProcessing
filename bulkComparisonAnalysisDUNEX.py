""""this scrips assumse a load from """
import pickle
import numpy as np
from matplotlib import pyplot as plt
import os
from testbedutils import sblib as sb
from matplotlib.colors import LogNorm
fpath = "plots/DunexPres"
# load the pickle with my data
logStats = pickle.load(open("logStats.pkl", 'rb'))
# unpack my data
datestring, residualsE, surveyE, newXE, crawlerE, profileNumber, stats, newX, survey, crawler, residuals = [], [], [], \
                                        [], [], [], [],  [], [], [], []
for i in range(len(logStats)):
    datestring.append(logStats[i][0])
    profileNumber.append(logStats[i][1])
    stats.append(logStats[i][2])
    newX.append(logStats[i][3])
    survey.append(logStats[i][4])
    crawler.append(logStats[i][5])
    surveyE.extend(logStats[i][4])
    crawlerE.extend(logStats[i][5])
    print(f" date {datestring[i]}, profileNumber {profileNumber[i]}")
    residuals.append(crawler[i] - survey[i])
    residualsE.extend(crawler[i] - survey[i])
    newXE.extend(logStats[i][3])

#calculate overall statistis
statsT = sb.statsBryant(surveyE, crawlerE)

# binned residuals and 1:1
plt.figure();
ax1 = plt.subplot2grid((2,1), (0,0))
ax1.plot(surveyE, crawlerE, '.')
ax1.plot([-3, 3], [-3, 3], 'k-')
ax1.set_xlabel('survey elevation [m]')
ax1.set_ylabel('crawler elevation [m]')
ax2 = plt.subplot2grid((2,1), (1,0))
ax2.hist(residualsE, bins=100);
ax2.set_xlabel('residuals [m]')
plt.savefig(os.path.join(fpath, "binnedResiduals"))
plt.close()

#exploreable plots
plt.figure();
plt.hist2d(surveyE, crawlerE, bins=200, cmap='YlOrBr', norm=LogNorm())
a=plt.colorbar();
a.set_label('counts');
plt.title('depth comparison')
plt.plot([-3, -3], [3,3], 'k-')
plt.xlabel('survey elevation NAVD88 [m]')
plt.ylabel('crawler elevation NAVD88 [m]')
plt.savefig(os.path.join(fpath, "SkillComparison"))
plt.close()

plt.figure();
plt.hist2d(surveyE, residualsE, cmap='YlOrBr',  bins=100, norm=LogNorm());
plt.plot([-3,3], [0,0], 'k-');
plt.xlabel('depth[m]');
plt.ylabel('residuals');
a=plt.colorbar();
a.set_label('counts')
plt.title('residuals of depth')
plt.savefig(os.path.join(fpath, "depthResiduals"))
plt.close()


plt.figure()
plt.hist2d(newXE, residualsE, bins=100, cmap='YlOrBr', norm=LogNorm())
a = plt.colorbar()
a.set_label('counts')
plt.plot([50, 275], [0,0], 'k-')
plt.xlabel('xFRF')
plt.ylabel('residuals')
plt.title('residuals of cross-shore')
plt.savefig(os.path.join(fpath, "crossShoreResiduals"))
plt.close()

# 1:1 plot


# heat map 1:1
plt.savefig(os.path.join(fpath, "binnedResiduals"))
plt.close()
# bin residuals in cross-shore location

