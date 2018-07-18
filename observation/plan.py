#######################################################################
# Ethasham Ahmed, 2018
# King's College London
#
# This is a program for observation planning using astropy.
# The program sorts through a given catalog of astronomical bodies and finds
# best targets for observation based on altitude and airmass.
#######################################################################
#importing necessary modules
from __future__ import print_function
import numpy as np
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_constellation, ICRS, Angle, Latitude, Longitude
from astropy.time import Time
from urllib.parse import urlencode
from urllib.request import urlretrieve
from astropy.table import Table
from IPython.display import Image
import os

# Sets observation time from user input
t = input('Enter time in (yyyy-m-d hh:mm:ss) format: ')
time = Time(t, scale='utc')
# Sets observation location
strand = EarthLocation(lat=51.511136*u.deg,lon=-0.115796*u.deg,height=20*u.m)

# This is the main function of the code which runs all the other functions
def main():
    q, data = observation()
    output(data)
    print('The best objects for viewing: ')
    print('Name\t\tAltitude\t\tAzimuth\t\tAirmass\t\tconstellation')
    # Data for the best 15 targets are shown, if there arent more than 15
    # targets available to observe then code prints however many number of targets that are available.
    if len(data) > 15:
        for n in range(15):
            print('{0}\t{1}\t{2}\t{3}\t{4}'.format(data[n][0], data[n][1], data[n][2],data[n][3], data[n][4]))
        if input('Do you want to download images of the objects? (y/n) ') == 'y':
            # Creates a new folder to store images
            if not os.path.exists("images"):
                os.makedirs("images")
            print('Downloading')
            for n in range(15):
                get_image(q[n],n)

    else:
        for n in range(len(data)):
            print('{0}\t{1}\t{2}\t{3}\t{4}'.format(data[n][0], data[n][1], data[n][2],data[n][3], data[n][4]))
        if input('Do you want to download images of the objects? (y/n) ') == 'y':
            # Creates a new folder to store images
            if not os.path.exists("images"):
                os.makedirs("images")
            print('Downloading')
            for n in range(15):
                get_image(q[n],n)

# This function gathers all the necessary data from a given catalog
def observation():
    # Limits the number of the object used from the catalog
    Vizier.ROW_LIMIT = 100 # setting the value to -1 will remove limit.
    ans = input('Do you want to search for a catalog or an object? ')
    # Uses Vizier to search for catalogs of astronomical bodies.
    if ans == 'catalog':
        catalogs = Vizier.get_catalogs(input('Enter catalog name: '))
    elif ans == 'object':
        catalogs = Vizier.query_object(input('Enter object name: '))

    print(catalogs)
    # Stores the catalog required from the search results.
    result = catalogs[int(input('Enter catalog number: '))]
    # Prints format of the catalogs data
    print(result.info)
    print(result)
    r = int(input('Enter column number for Right Ascension(Ra): '))
    d = int(input('Enter coloumn number for Declination(Dec): '))
    # Stores coordinates in the ICRS frame.
    q = [SkyCoord(result[n][r],result[n][d], 'icrs', unit="deg") for n in range(len(result))]
    data = []
    for i in range(len(q)):
        name = result[i][0]
        # Transforms the coordinates to Altitude and Azimuth using observation time and location
        altaz = q[i].transform_to(AltAz(obstime=time, location=strand))
        # Retrives Constellation name using Altitude and Azimuth
        const = get_constellation(altaz)
        # Filters out objects that are too high
        if altaz.alt.degree <= 75:
            data.append([name, altaz.alt.degree, altaz.az.degree,altaz.secz, const])
    # Stores all necessary data obtained into an array sorted by altitude.
    data = sorted(data, key=lambda alt:alt[1], reverse=True)
    return q, data

# This function writes all the necessary data to a text file.
def output(data):
    out_file = open('targets.txt', 'w')
    out_file.write('Name\t\tAltitude\t\tAzimuth\t\t\tAirmass\t\tconstellation\n')
    for n in range(len(data)):
        out_file.write('%s	%s	%s	%s  %s\n' % (data[n][0], data[n][1], data[n][2],data[n][3], data[n][4]))
    out_file.close()

# Retrives images for the best targets.
def get_image(q,n):
    impix = 1024
    imsize = 12*u.arcmin
    # This link is for the SDSS database
    sourceurl = 'http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx'
    # Uses Ra and Dec to retrive the required image
    query_string = urlencode(dict(ra=q.ra.deg,
                                  dec=q.dec.deg,
                                  width=impix, height=impix,
                                  scale=imsize.to(u.arcsec).value/impix))
    url = sourceurl + '?' + query_string
    # Downloads the image
    image_name = 'images/Object'+str(n)+'.jpg'
    urlretrieve(url, image_name)
    Image(image_name)

main()
