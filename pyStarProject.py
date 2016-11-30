#!/usr/bin/env python
import numpy
import dodecahedron
import csv

class Star:
    """ class to store a star's data from the catalogue database """
    id = -1 # The database primary key.
    ra = 0 # The star's right ascension, for epoch and equinox 2000.0.
    dec = 0 # The star's declination, for epoch and equinox 2000.0.
    name = "" # A common name for the star,
    mag = 0 # The star's apparent visual magnitude.
    def __init__(self, db):
        self.id = db['id']
        self.ra = db['ra']
        self.dec = db['dec']
        self.name = db['proper']
        self.mag = db['mag']
        
        
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
            if starreader.line_num > 200:
                break
    

            
