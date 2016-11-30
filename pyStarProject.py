#!/usr/bin/env python
import numpy
import dodecahedron
import csv

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
    print("Dodecahedron with vertices on poles")
    vertexDodec = dodecahedron.Dodecahedron(withFacesOnPoles=False)
    for i in range(12):
        faceVec = vertexDodec.getFaceCtr(i)
        print("Face %2d: %s" % (i, faceVec))

    starlist = []
    with open('HYG-Database/hygdata_v3.csv', 'r') as csvfile:
        starreader = csv.DictReader(csvfile, delimiter=',')
        for row in starreader:
            starlist.append(Star(row))
            if starreader.line_num > 2000:
                break

    x = np.array([s.ra for s in starlist])
    plt.hist(x, 50)
    plt.show()
