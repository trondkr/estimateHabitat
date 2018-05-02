# estimateHabitat

<h1>estimateHabitat</h1>
<a href="https://badge.fury.io/gh/trondkr%2FestimateHabitat"><img src="https://badge.fury.io/gh/trondkr%2FestimateHabitat.svg" alt="GitHub version" height="18"></a>

<a href="https://codeclimate.com/github/trondkr/estimateHabitat/maintainability"><img src="https://api.codeclimate.com/v1/badges/302e97e08f2f553f9616/maintainability" /></a>

EstimateHabitat is a selection of files used to calculate the changes in light conditions in the upper part of the water column in the ocean accounting for thickness of sea-ice and snow andgorgraphical locations and time of the year. The files area a combination of Fortran and Python files compiled using Cython (setup.py). 

<h3>Latest updates</h3>
<ul>
<li>02.05.2018: Added the files to Github</li>
</ul>

<h3>Order of execution</h3>
<ul>
    <li>Compile: A cython function which is found in the file calculateLightUnderIce.pyx is
    compiled with: python setup.py build_ext --inplace</li>
    <li>Download the NorESM historical and future predictions files (any scenario) as input: {variable}_Omon_NorESM1-M_%s_r1i1p1_200601-210012_rectangular.nc and {variable}_OImon_NorESM1-M_historical_r1i1p1_185001-200512_rectangular.nc where variable is all of ["sit", "sic", "snd", "ialb", "tos"] </li>
    <li>calculateMaxLightEntireArctic.py is run first and creates an output file: Light_and_temperature_1850_2100_Arctic.nc</li>
    <li>Execute estimateHabitat.py to calculate the geographical changes in tos and light combined for decadal periods</li>
    <li></li>
</ul>
