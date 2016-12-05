#!/usr/bin/env python3
import csv                               ## for reading in our data files
import argparse
import logging                           ## for orderly print output

import math
import numpy as np
import dodecahedron

import matplotlib.pyplot as plt

class Star:
    """ class to store a star's data from the catalogue database """

    def __init__(self, idno, ra, dec, name, mag, highlight = False):
        self.id = int(idno) # The database primary key.
        self.ra = float(ra) # The star's right ascension, for epoch and equinox 2000.0.
        self.dec = float(dec)# The star's declination, for epoch and equinox 2000.0.
        self.name = name # A common name for the star
        self.mag = float(mag)  # The star's apparent visual magnitude.
        self.highlight = highlight
        # convert the ra and dec into angles in rad
        phi = (self.ra/24) * 2 * math.pi
        rho = (self.dec/90) * math.pi
        self.vec = np.array(np.zeros(3))
        self.vec[0] = math.sin(rho) * math.cos(phi) #x
        self.vec[1] = math.sin(rho) * math.sin(phi) #y
        self.vec[2] = math.cos(rho) #z

    @classmethod
    def fromcatalogue(cls, db):
        "Initialize Star from a catalogue db entry"        
        return cls(idno = db['id'], ra = db['ra'], dec = db['dec'],
                   name = db['proper'], mag = db['mag'])
    @classmethod
    def fromstar(cls, star):
        "Initialize Star from another star"        
        return cls(idno = star.id, ra = star.ra, dec = star.dec,
                   name = star.name, mag = star.mag, highlight = star.highlight)

class PinnedStar(Star):
    """ Class for pinning a star's coordinates onto dodecahedron """
    
    def __init__(self, star, dodecah):
        Star.__init__(self, idno = star.id,
                        ra = star.ra, dec = star.dec,
                        name = star.name, mag = star.mag, highlight = star.highlight)
        # calculate which face to project on
        self.faceid = dodecah.getFaceInd(self.vec)
        # determine cartesian coordinates for star on dodecahedron
        facectr = dodecah.getFaceCtr(self.faceid)
        self.dodecvec = self.vec * (1/np.dot(self.vec, facectr))
        
if __name__ == "__main__":
    ## set up some print-out routines (logging)
    FORMAT = '%(asctime)s %(name)s:line %(lineno)-4d %(levelname)-8s %(message)s'
    logging.basicConfig(format=FORMAT)
    log = logging.getLogger('pyStarProject') ## set up logging
    log.setLevel("DEBUG")
    
    # argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--nentries", help="max number of entries to read from the catalogue",
                        type=int, default=-1, action="store")
    parser.add_argument("--mag", help="largest apparent magnitude to still consider for selecting stars from catalogue",
                        type=float, default=6.5, action="store")
    parser.add_argument("--highlight", dest='highlight', help='what constellations (given by their abbreviation) to hightlight', nargs='+')
    args = parser.parse_args()

    # read in star catalogue data
    starlist = []
    with open('HYG-Database/hygdata_v3.csv', 'r') as csvfile:
        starreader = csv.DictReader(csvfile, delimiter=',')
        for row in starreader:
            s = Star.fromcatalogue(row)
            if s.mag < args.mag and s.mag > -26: # filter low-visibility stars and sun
                if args.highlight and row['con'] in args.highlight:
                    log.debug("Found constellation '{}' to highlight".format(row['con']))
                    s.highlight = True
                if args.highlight and s.name in args.highlight:
                    log.debug("Found star '{}' to highlight".format(s.name))
                    s.highlight = True
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
    maglist = np.array([s.mag for s in starlist])
    plt.hist(maglist, 50)
    plt.xlabel("mag")

    # marker size for 3D plot should depend on magnitude
    magmin = np.min(maglist)
    magmax = np.max(maglist)
    msizemax = 100
    msizemin = 1
    def msize(mag):
        return msizemax + (mag - magmin)*((msizemin - msizemax)/(magmax - magmin))

    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter([s.vec[0] for s in starlist],[s.vec[1] for s in starlist], [s.vec[2] for s in starlist],
               s=[msize(s.mag) for s in starlist], color=["red" if s.highlight else "blue" for s in starlist])

    # create dodecahedron
    print("Dodecahedron with vertices on poles")
    Dodec = dodecahedron.Dodecahedron(withFacesOnPoles=False)
    for i in range(12):
        faceVec = Dodec.getFaceCtr(i)
        print("Face %2d: %s, length^2 = %f" % (i, faceVec, np.dot(faceVec, faceVec)))
        
    # create list of stars pinned on dodec's faces
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    pinnedstarlist = [PinnedStar(s, Dodec) for s in starlist]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # draw dodecahedron with pinned stars
    verts = []
    for i in range(12):
        faceVec = Dodec.getFaceCtr(i)
        thisverts = Dodec.getVertices(i)
        norm = 1/np.dot(thisverts[0], faceVec)
        for i in thisverts:
            i *= 0.95*norm # normalize vertices for displaying dodecahedron with right size
        verts.append(thisverts)
        ax.scatter(faceVec[0], faceVec[1], faceVec[2], c='r')
    ax.add_collection3d(Poly3DCollection(verts))
    ax.scatter([s.dodecvec[0] for s in pinnedstarlist],[s.dodecvec[1] for s in pinnedstarlist],
               [s.dodecvec[2] for s in pinnedstarlist],
               s=[msize(s.mag) for s in pinnedstarlist],
               color=["red" if s.highlight else "blue" for s in pinnedstarlist])
    
    # for one face, construct a 2D image of the face on the face's local xy-plane
    # orientation: local x in global z direction
    for i in range(12):
        # construct unit vectors in the face's plane
        vecface = Dodec.getFaceCtr(i)
        veclocy = np.cross(vecface ,np.array([0,0,1]))
        veclocy *= 1/math.sqrt(np.dot(veclocy,veclocy))
        veclocx = np.cross(vecface,veclocy)
        veclocx *= 1/math.sqrt(np.dot(veclocx, veclocx))
        # verify vector orientation
        if np.dot(veclocx, [0,0,1]) < 0:
            print ("Face {} x vector reorientation needed".format(i))
            veclocx *= -1
        if np.dot(np.cross(veclocx, veclocy), vecface) < 0:
            print ("Face {} y vector reorientation needed".format(i))
            veclocy *= -1
            
        facestarlist = [s for s in pinnedstarlist if s.faceid == i]
        plt.figure()
        plt.scatter([np.dot(s.vec-vecface,veclocx) for s in facestarlist],
                    [np.dot(s.vec-vecface,veclocy) for s in facestarlist],
                    s=[msize(s.mag) for s in facestarlist],
                    color=["red" if s.highlight else "blue" for s in starlist])
        plt.axis((-1,1,-1,1))
    
    ## final commands to show resulting plots (not needed in an interactive session)
    plt.draw() ## <- shows the plots
    plt.pause(1) # <-------
    input("<Hit Enter To Close>")
    plt.close('all')

