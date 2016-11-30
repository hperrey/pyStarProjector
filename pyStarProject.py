#!/usr/bin/env python
import numpy
import dodecahedron
import csv

if __name__ == "__main__":
    print("Dodecahedron with vertices on poles")
    vertexDodec = dodecahedron.Dodecahedron(withFacesOnPoles=False)
    for i in range(12):
        faceVec = vertexDodec.getFaceCtr(i)
        print("Face %2d: %s" % (i, faceVec))

    with open('HYG-Database/hygdata_v3.csv', 'r') as csvfile:
        starreader = csv.DictReader(csvfile, delimiter=',')
        for row in starreader:
            print (row)
            if starreader.line_num > 20:
                break
