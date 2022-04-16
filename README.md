# rmt-sns : software for remote sensing tasks

## Build

Build and test on ubuntu 18; first install dependencies:

> sudo ./setup.sh

then simply make in this directory, which should create three apps in the bin dir. 

## Examples

To create images from apt wav files:

> bin/aptdec data/NOAA1520200811-143715.wav

You can model the camera view for this satellite image, including country coordinates and lon/lat using the 'ig' app. For example, eye just above Queensland:

> bin/ig --dr 0.5  --wp data/countries.kml --cov data/test_image2.png --d 10000 --dp 2 &> /dev/null

which creates output in the data directory, showing country outlines and lon/lat grid points. Run app with --help to see all options.

This application works primarily with the wm (web mercator) application. In the following execution, a model for the image grid is created (maps pixel location to lon/lat), which is piped to 'wm', then used to create tiles from an input image:

>  bin/ig --dr 0.04  --wp data/countries.kml --dp 2 2> /dev/null | bin/wm --in data/test_image2.png --z 6 --out ./data/

Output tiles can easily be qc'd against open street maps eg., data/6_57_34.png <=> https://tile.openstreetmap.org/6/57/34.png

WJB 08/21
   
