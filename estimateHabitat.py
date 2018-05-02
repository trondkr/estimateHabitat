from pylab import *
import os, datetime
import numpy as np
from netCDF4 import Dataset, date2num, num2date
import mpl_util
import sys
import pandas as pd
from operator import and_

from mpl_toolkits.basemap import Basemap, interp, shiftgrid, addcyclic
import brewer2mpl
import calendar
from scipy.stats.mstats import gmean

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime.datetime(2017, 12, 20)
__modified__ = datetime.datetime(2017, 12, 20)
__version__ = "1.1"
__status__ = "Development, 20.12.2017"

"""This script calculates change in habitat from historical values. Habitat is defined
as a range of temperature and light in the water column for 4 seasons. 

This script requires the output from running:
calculateMaxLightEntireArctic.py

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

    # turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')

    # now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()
    if right:
        ax.yaxis.tick_right()


def plotTimeseries(ts, myvar, season):
    ts_annual = ts.resample("A")
    ts_quarterly = ts.resample("Q")
    ts_monthly = ts.resample("M")

    # Write data to file
    mypath = "%s_annualaverages.csv" % (myvar)
    if os.path.exists(mypath): os.remove(mypath)
    ts.to_csv(mypath)
    print(("Wrote timeseries to file: %s" % (mypath)))

    red_purple = brewer2mpl.get_map('RdPu', 'Sequential', 9).mpl_colormap
    colors = red_purple(np.linspace(0, 1, 12))
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    # for mymonth in xrange(12):
    #    ts[(ts.index.month == mymonth + 1)].plot(marker='o', color=colors[mymonth], markersize=5, linewidth=0,
    #                                            alpha=0.8)
    # ts_annual.plot(marker='o', color="#FA9D04", linewidth=0, alpha=1.0, markersize=7, label="Annual")
    remove_border(top=False, right=False, left=True, bottom=True)
    ts.resample("M").mean().plot(style="r", marker='o', linewidth=1, label="Monthly")
    ts.resample("A").mean().plot(style="b", marker='o', linewidth=2, label="Annual")

    # legend(loc='best')
    if myvar == "light":
        ylabel(r'Light (W m$^{-2})$')

    if myvar == "temp":
        ylabel(r'Temp ($^{o}$C)')

    plotfile = 'figures/timeseries_' + str(season) + '_' + str(myvar) + '.png'
    plt.savefig(plotfile, dpi=300, bbox_inches="tight", pad_inches=0)

    plt.show()


def getData(infile):
    if os.path.exists(infile):
        try:
            cdf = Dataset(infile)
            print(("Opened inputfile: %s" % (infile)))
        except:
            print(("Unable to  open file: %s" % (infile)))
            sys.exit()

        temp = cdf.variables["tos"][:]
        light = cdf.variables["light"][:]
        times = cdf.variables["time"][:]
        longitude = cdf.variables["longitude"][:]
        latitude = cdf.variables["latitude"][:]

        dates = num2date(times, "days since 1948-01-01 00:00:00", calendar="365_day")
        print("Extracted time-steps starting in %s and ending in %s" % (dates[0], dates[-1]))

        return temp, light, dates, longitude, latitude


def getStartAndEndIndex(startYear, endYear, dates):
    startIndex = -9;
    endIndex = -9
    for dateIndex, JD in enumerate(dates):

        if JD.year == startYear:
            startIndex = dateIndex
        if JD.year == endYear:
            endIndex = dateIndex
    if startIndex == -9 or endIndex == -9:
        print(("Unable to find indexes for start %s and end %s years", startYear, endYear))
        sys.exit()

    print(("=> Period %s to %s" % (dates[startIndex], dates[endIndex])))
    return startIndex, endIndex


def createDecadalAverages(seasonsArray, seasons, latitude, longitude, dates, periods):
    # New array to store values will contain seasonal values of tos and light as well as standard deviations
    # within each period of each.
    decadalArray = np.zeros((6, len(seasons), len(periods)-1, len(latitude), len(longitude)))


    for seasonIndex, season in enumerate(seasons):

        for periodIndex in range(len(periods)-1):
            startIndex, endIndex = getStartAndEndIndex(periods[periodIndex], periods[periodIndex + 1], dates)
            # 0= TOS, 1=STD(TOS), 2 = LIGHT, 3=STD(LIGHT), 4=Normalized tos, 5=normalized light

            temp_period = np.squeeze(seasonsArray[0, seasonIndex, startIndex:endIndex, :, :])
            light_period = np.squeeze(seasonsArray[1, seasonIndex, startIndex:endIndex, :, :])
            norm_temp_period = np.squeeze(seasonsArray[2, seasonIndex, startIndex:endIndex, :, :])
            norm_light_period = np.squeeze(seasonsArray[3, seasonIndex, startIndex:endIndex, :, :])

            decadalArray[0, seasonIndex, periodIndex, :, :] = np.squeeze(
                np.ma.mean(temp_period, axis=0))
            decadalArray[1, seasonIndex, periodIndex, :, :] = np.squeeze(
                np.ma.std(temp_period, axis=0))
            decadalArray[2, seasonIndex, periodIndex, :, :] = np.squeeze(
                np.ma.mean(light_period, axis=0))
            decadalArray[3, seasonIndex, periodIndex, :, :] = np.squeeze(
                np.ma.std(light_period, axis=0))
            decadalArray[4, seasonIndex, periodIndex, :, :] = np.squeeze(
                np.ma.mean(norm_temp_period, axis=0))
            decadalArray[5, seasonIndex, periodIndex, :, :] = np.squeeze(
                np.ma.mean(norm_light_period, axis=0))

            decadalArray = ma.masked_invalid(decadalArray)
            print(("MEAN temp %s" % (np.ma.mean(decadalArray[0, seasonIndex, periodIndex, :, :]))))
            print(("MEAN light %s" % (np.ma.mean(decadalArray[1, seasonIndex, periodIndex, :, :]))))
            print(("MEAN std temp %s" % (np.ma.mean(decadalArray[2, seasonIndex, periodIndex, :, :]))))
            print(("MEAN std light %s" % (np.ma.mean(decadalArray[3, seasonIndex, periodIndex, :, :]))))
            print(("MEAN norm temp %s" % (np.ma.mean(decadalArray[4, seasonIndex, periodIndex, :, :]))))
            print(("MEAN norm light %s" % (np.ma.mean(decadalArray[5, seasonIndex, periodIndex, :, :]))))

        #  test=np.rot90(np.flipud(np.squeeze(decadalArray[4, seasonIndex, periodIndex, :, :])), 3)
        #  plotMap(longitude, latitude, test, periods[periodIndex], "test", seasons[seasonIndex])

        # test = np.rot90(np.flipud(np.squeeze(decadalArray[2, seasonIndex, periodIndex, :, :])), 3)
        # plotMap(longitude, latitude, test, periods[periodIndex], seasons[seasonIndex])

    return decadalArray, periods


def printStatistics(decadalArray, seasons, periods):
    for seasonIndex, season in enumerate(seasons):
        for periodIndex in range(len(periods) - 1):

            if periodIndex > 0:
                avgPeriod = np.mean(decadalArray[0, seasonIndex, periodIndex, :, :])
                avgClimPeriod = np.mean(decadalArray[0, seasonIndex, 0, :, :])
                change = ((avgPeriod - avgClimPeriod) / avgClimPeriod) * 100.


def estimateChangedEcosystem(decadalArray, seasons, periods):
    print("seasons {}".format(seasons))
    print("periods {}".format(periods))

    # Size: [variables, seasons, period diffs, lat, long]
    estimateChangeArray = np.zeros((3, np.shape(decadalArray[0, :, 0, 0, 0])[0],
                                    np.shape(decadalArray[0, 0, :, 0, 0])[0]-1,
                                    np.shape(decadalArray[0, 0, 0, :, 0])[0],
                                    np.shape(decadalArray[0, 0, 0, 0, :])[0]))
    print("decadalArray {}".format(np.shape(decadalArray)))


    for seasonIndex, season in enumerate(seasons):
        for periodIndex in range(1, np.shape(decadalArray[0, 0, :, 0, 0])[0]):
            for lat in range(np.shape(decadalArray[0, 0, 0, :, 0])[0]):
                for lon in range(np.shape(decadalArray[0, 0, 0, 0, :])[0]):

                    tos = decadalArray[4, seasonIndex, periodIndex, lat, lon]
                    light = decadalArray[5, seasonIndex, periodIndex, lat, lon]

                    tosclim = decadalArray[4, seasonIndex, 0, lat, lon]
                    lightclim = decadalArray[5, seasonIndex, 0, lat, lon]

                    QCTOS = 0
                    QCLIGHT = 0
                    print("TEST: {} - tos {} light {} toclim {} lightclim {} periods: {}".format(periods[periodIndex], tos,light, tosclim, lightclim, np.shape(decadalArray[0, :, 0, 0, 0])[0]))

                    stds = np.arange(0.1, 2.0, 0.1)  # [0.5, 1.0, 1.5, 2.0, 3.0]
                    for std in stds:
                        if abs(tos - tosclim) > std:
                            QCTOS = std
                        # print("TOS: {}".formst(tos, tosclim, )
                        if abs(light - lightclim) > std:
                            QCLIGHT = std
                    print("form: {} sea: {} per: {}".format(np.shape(estimateChangeArray), seasonIndex, periodIndex))
                    estimateChangeArray[0, seasonIndex, periodIndex-1, lat, lon] = QCTOS
                    estimateChangeArray[1, seasonIndex, periodIndex-1, lat, lon] = QCLIGHT
                    # ß  print(QCTOS, QCLIGHT, tos, tosclim)
                    estimateChangeArray[2, seasonIndex, periodIndex-1, lat, lon] = gmean([QCTOS, QCLIGHT])

    estimateChangeArray = np.ma.masked_invalid(estimateChangeArray)

    return estimateChangeArray


def plotMap(lon, lat, mydata, period, qctype, season):
    plt.figure(figsize=(12, 12), frameon=False)
    mymap = Basemap(projection='npstere', lon_0=0, boundinglat=50)

    llat, llon = np.meshgrid(lat, lon)
    print("Plotting season: {} for period {}".format(season, period))
    x, y = mymap(llon, llat)
    print(np.min(mydata), np.max(mydata))
    levels = np.arange(np.min(mydata), np.max(mydata), 0.01)
    # levels = [0, 1]

    CS1 = mymap.contourf(x, y, mydata, levels,
                         cmap=mpl_util.LevelColormap(levels, cmap=cm.RdBu_r),
                         extend='max')

    mymap.drawparallels(np.arange(-90., 120., 15.), labels=[1, 0, 0, 0])  # draw parallels
    mymap.drawmeridians(np.arange(0., 420., 30.), labels=[0, 1, 0, 1])  # draw meridians

    mymap.drawcoastlines()
    mymap.drawcountries()
    mymap.fillcontinents(color='grey', alpha=0.2)
    plt.colorbar(CS1, shrink=0.5)
    title('QC:' + str(period) + ' season:' + str(season))

    CS1.axis = 'tight'
    if not os.path.exists("Figures"):
        os.mkdir("Figures/")
    plotfile = 'figures/map_qc_' + str(period) + '_season_' + str(season) + '.png'
    # plt.show()

    plt.savefig(plotfile, dpi=100, bbox_inches='tight')
    plt.clf()
    plt.close()


def createSeasonArrays(temp, light, dates, longitude, latitude):
    periods = [1950, 2000, 2050, 2100]

    winter = ["Jan", "Feb", "Mar"]
    spring = ["Apr", "May", "Jun"]
    summer = ["Jul", "Aug", "Sep"]
    autumn = ["Oct", "Nov", "Dec"]

    seasons = [winter, spring, summer, autumn]
    seasonNames = ["winter", "spring", "summer", "autumn"]

    count = 0
    last_year = -9
    years = []
    for d in dates:
        if d.year > last_year and d.month == 1:
            last_year = d.year

            years.append(datetime.datetime(d.year, d.month, d.day))
            count += 1
    count += 1
    print(("Timeseries contains %s years" % (count)))

    seasonsArray = np.zeros((4, len(seasons), count, len(latitude), len(longitude)))

    for seasonIndex, season in enumerate(seasons):
        tindex = 0
        seasonName = seasonNames[seasonIndex]
        tempT = np.zeros((len(latitude), len(longitude)))
        tempL = np.zeros((len(latitude), len(longitude)))
        yearFinished = -9
        counter = 0

        for dateIndex, JD in enumerate(dates):

            if calendar.month_abbr[JD.month] in season:
                tempT = tempT + temp[dateIndex, :, :]
                tempL = tempL + light[dateIndex, :, :]
                counter += 1
            else:
                if yearFinished != JD.year:
                    if counter > 0:
                        seasonsArray[0, seasonIndex, tindex, :, :] = tempT / counter * 1.0
                        seasonsArray[1, seasonIndex, tindex, :, :] = tempL / counter * 1.0

                    tindex += 1
                    tempT = tempT * 0.0
                    tempL = tempL * 0.0
                    counter = 0
                    yearFinished = JD.year

        # Create a set of normalized tos and light data
        for j in range(len(latitude)):
            for i in range(len(longitude)):
                tosN = seasonsArray[0, seasonIndex, :, j, i]
                lightN = seasonsArray[1, seasonIndex, :, j, i]
                seasonsArray[2, seasonIndex, :, j, i] = (tosN - np.ma.min(tosN)) / (np.ma.max(tosN) - np.ma.min(tosN))
                seasonsArray[3, seasonIndex, :, j, i] = (lightN - np.ma.min(lightN)) / (
                        np.ma.max(lightN) - np.ma.min(lightN))

        createTimeseriesPlot = False
        if createTimeseriesPlot:
            ll = []
            tt = []

            for i in range(len(years)):
                tt.append(np.ma.mean(seasonsArray[0, seasonIndex, i, :, :]))
                ll.append(np.ma.mean(seasonsArray[1, seasonIndex, i, :, :]))

            tsl = pd.Series(ll, years)
            plotTimeseries(tsl, "light", seasonName)

            tst = pd.Series(tt, years)
            plotTimeseries(tst, "temp", seasonName)

    decadalArray, periods = createDecadalAverages(seasonsArray, seasons, latitude, longitude, years, periods)
    estimateChangeArray = estimateChangedEcosystem(decadalArray, seasons, periods)

    plotvars = ["tosQC", "lightQC", "combinedQC"]

    for seasonIndex, season in enumerate(seasons):

        # If you have 4 reference dates 1950, 2000, 2050, 2100, then the periodIndex here refers to
        # the difference periods between each date. E.g. 0 = 1950-2000, 1=2000-2050 etc.
        # One less than perdiods definition.
        for periodIndex in range(np.shape(estimateChangeArray[0,0,:,0,0])[0]):
            for i in range(len(plotvars)):
                mydata = np.squeeze(estimateChangeArray[i, seasonIndex, periodIndex, :, :])

                mydata = np.rot90(np.flipud(mydata), 3)

                periodName = "{}-{}-{}".format(plotvars[i], periods[periodIndex], periods[periodIndex + 1])
                print("Periodname {}".format(periodName))
                plotMap(longitude, latitude, mydata, periodName, plotvars[i], seasons[seasonIndex])


print(('Python %s on %s' % (sys.version, sys.platform)))

infile = "Light_and_temperature_1850_2100_Arctic.nc"
temp, light, dates, longitude, latitude = getData(infile)

createSeasonArrays(temp, light, dates, longitude, latitude)
