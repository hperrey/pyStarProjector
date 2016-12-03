#!/usr/bin/env python
import numpy
import dodecahedron
import csv
import argparse

import math
import numpy as np

import matplotlib.pyplot as plt

class Star:
    """ class to store a star's data from the catalogue database """
    id = -1 # The database primary key.
    ra = 0 # The star's right ascension, for epoch and equinox 2000.0.
    dec = 0 # The star's declination, for epoch and equinox 2000.0.
    name = "" # A common name for the star,
    mag = 0 # The star's apparent visual magnitude.
    def __init__(self, db):
        self.id = int(db['id'])
        self.ra = float(db['ra'])
        self.dec = float(db['dec'])
        self.name = db['proper']
        self.mag = float(db['mag'])
        
        
if __name__ == "__main__":
    # argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--nentries", help="max number of entries to read from the catalogue",
                        type=int, default=-1, action="store")
    parser.add_argument("--mag", help="largest apparent magnitude to still consider for selecting stars from catalogue",
                        type=int, default=6.5, action="store")
    args = parser.parse_args()

    # dodecahedron test
    print("Dodecahedron with vertices on poles")
    vertexDodec = dodecahedron.Dodecahedron(withFacesOnPoles=False)
    for i in range(12):
        faceVec = vertexDodec.getFaceCtr(i)
        print("Face %2d: %s" % (i, faceVec))

    # read in star catalogue data
    starlist = []
    with open('HYG-Database/hygdata_v3.csv', 'r') as csvfile:
        starreader = csv.DictReader(csvfile, delimiter=',')
        for row in starreader:
            s = Star(row)
            if s.mag < args.mag:
                starlist.append(s)
            if args.nentries > 0 and len(starlist) > args.nentries:
                print ("Reached requested maximum number of events: {}".format(args.nentries))
                break

    print ("Loaded {} stars into memory".format(len(starlist)))
            
    # plot some of the relevant star data
    plt.figure()
    plt.hist(np.array([s.dec for s in starlist]), 50)
    plt.xlabel("Dec")
    plt.figure()
    plt.hist(np.array([s.ra for s in starlist]), 50)
    plt.xlabel("RA")
    plt.figure()
    plt.hist(np.array([s.mag for s in starlist]), 50)
    plt.xlabel("mag")
    
    locations = []
    for star in starlist:
        # convert the ra and dec into angles in rad
        phi = (star.ra/24) * 2 * math.pi
        rho = (star.dec/90) * math.pi
        vec = np.array([math.sin(rho) * math.cos(phi), #x
                       math.sin(rho) * math.sin(phi), #y
                       math.cos(rho)] #z
                       )
        locations.append(vec)

    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter([l[0] for l in locations],[l[1] for l in locations], [l[2] for l in locations])
    plt.show()
