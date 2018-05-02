from pylab import *
import matplotlib.dates as mdates

import os, sys, datetime, string
import numpy as np
from netCDF4 import Dataset
import time
import numpy.ma as ma
import pandas as pd
import brewer2mpl

from mpl_toolkits.basemap import Basemap, interp, shiftgrid, addcyclic

__author__   = 'Trond Kristiansen'
__email__    = 'trond.kristiansen@imr.no'
__created__  = datetime.datetime(2010, 1, 16)
__modified__ = datetime.datetime(2014, 3, 4)
__version__  = "1.1"
__status__   = "Development, 16.01.2010, 14.04.2010, 04.03.2014"

"""This script plots the light averaged per year from running calculateMaxLight.py
"""



def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
    """
    Minimize chartjunk by stripping out unnecesasry plot borders and axis ticks
    
    The top/right/left/bottom keywords toggle whether the corresponding plot border is drawn
    """
    ax = axes or plt.gca()
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)

    #turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')

    #now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()
    if right:
        ax.yaxis.tick_right()

def printStats(mydata):

    for d in mydata:
        print(mydata)
def plotTimeseries(ts,myvar):

    ts_annual = ts.resample("A")
    ts_10A = ts.resample("10A")
    ts_20A = ts.resample("20A")
#    printStats(ts_10A)

    # Write data to file
    mypath="%s_annualaverages.csv"%(myvar)
    if os.path.exists(mypath):os.remove(mypath)
    ts.to_csv(mypath)
    print("Wrote timeseries to file: %s"%(mypath))

    red_purple = brewer2mpl.get_map('RdPu', 'Sequential', 9).mpl_colormap
    colors = red_purple(np.linspace(0, 1, 12))
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
   # for mymonth in xrange(12):
        #ts[(ts.index.month==mymonth+1)].plot(marker='o', color=colors[mymonth],markersize=5, linewidth=0,alpha=0.8)
    #hold(True)
    ts_annual.plot(color="#FA9D04", linewidth=0.5,alpha=1.0, label="")
   
    ts_annual.plot(marker='o', color="#FA9D04", linewidth=0,alpha=1.0, markersize=5, label="Annual")
    ts_10A.plot(color="r", linewidth=0.5,alpha=1.0, label="10 years")
    ts_20A.plot(color="m", linewidth=0.5,alpha=1.0, label="20 years")
   
   # remove_border(top=False, right=False, left=True, bottom=True)
    #ts_monthly.plot(style="r", marker='o', linewidth=1,label="Monthly")

    # legend(loc='best')
    ylabel(r'Light (W m$^{-2})$')

    plotfile='figures/timeseries_'+str(myvar)+'.png'
    plt.savefig(plotfile,dpi=300,bbox_inches="tight",pad_inches=0)
    print('Saved figure file %s\n'%(plotfile))
    plt.show()


infile="light_annualaverages.csv"

f=open(infile,'r')
lines=f.readlines()
mydates=[]; myvalues=[]
counter=0
for line in lines:
    if counter>0:
        l=string.split(line,',')
        mydate=datetime.datetime.strptime(l[0], '%Y-%m-%d %H:%M:%S')
        mydates.append(mydate)
        myvalue=float(l[1].strip())
        myvalues.append(myvalue)
    
    counter+=1
f.close()
ts=pd.Series(myvalues,mydates)
print(ts)
plotTimeseries(ts,"light")
