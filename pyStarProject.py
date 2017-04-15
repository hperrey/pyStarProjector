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

    def __init__(self, idno, ra, dec, name, mag, vec = np.array(np.zeros(3)), highlight = False):
        self.id = int(idno) # The database primary key.
        self.ra = float(ra) # The star's right ascension, for epoch and equinox 2000.0.
        self.dec = float(dec)# The star's declination, for epoch and equinox 2000.0.
        self.name = name # A common name for the star
        self.mag = float(mag)  # The star's apparent visual magnitude.
        self.highlight = highlight
        if not vec.any():
            self.vec = _calccoordinates(self.ra, self.dec)
        else:
            self.vec = vec
        
    def shift(self, delta_ra, delta_dec):
        # shift star by delta RA and DEC
        self.dec += delta_dec
        self.ra += delta_ra
        # check that result is within scope of parameters
        self.dec = ((self.dec+90)%180)-90
        self.ra = self.ra%24
        self.vec = _calccoordinates(self.ra, self.dec)

    @classmethod
    def fromcatalogue(cls, db):
        "Initialize Star from a catalogue db entry"        
        return cls(idno = db['id'], ra = db['ra'], dec = db['dec'],
                   name = db['proper'], mag = db['mag'])
    @classmethod
    def fromstar(cls, star):
        "Initialize Star from another star"        
        return cls(idno = star.id, ra = star.ra, dec = star.dec,
                   name = star.name, mag = star.mag, vec = star.vec, highlight = star.highlight)
                   
    @classmethod
    def empty(cls):
        "Initialize Star with dummy info"        
        return cls(idno = -1, ra = 0, dec = 0,
                   name = "empty", mag = -20, highlight = True)

class PinnedStar(Star):
    """ Class for pinning a star's coordinates onto dodecahedron """
    
    def __init__(self, star, dodecah):
        Star.__init__(self, idno = star.id,
                        ra = star.ra, dec = star.dec,
                        name = star.name, mag = star.mag, 
                        vec = star.vec, highlight = star.highlight)
        # calculate which face to project on
        self.faceid = dodecah.getFaceInd(self.vec)
        # determine cartesian coordinates for star on dodecahedron
        facectr = dodecah.getFaceCtr(self.faceid)
        self.dodecvec = self.vec * (1/np.dot(self.vec, facectr))
        
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    # code from http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def _calccoordinates(ra, dec):
    """
     convert the ra and dec into angles (in rad) and return cartesian vector
     """
    phi = (ra/24) * 2 * math.pi
    rho = (dec/90) * math.pi
    vec = np.array(np.zeros(3))
    vec[0] = math.sin(rho) * math.cos(phi) #x
    vec[1] = math.sin(rho) * math.sin(phi) #y
    vec[2] = math.cos(rho) #z
    return vec
                     
def _findFacesSharingVertices(vertices, dodecahedron):
    '''returns which faces of a dodecahedron share the two given vertices'''
    faces = []
    for f in range(12):
        matches = 0
        for v in dodecahedron.getVertices(f):
            for ref in vertices:
                if np.equal(v, ref).all():
                    matches += 1
        if matches == len(vertices):
            faces.append(f)
    return faces

def main():
    ## set up some print-out routines (logging)
    FORMAT = '%(asctime)s %(name)s:line %(lineno)-4d %(levelname)-8s %(message)s'
    logging.basicConfig(format=FORMAT)
    log = logging.getLogger('pyStarProject') ## set up logging
    log.setLevel("DEBUG")
    
    # argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--nentries", help="max number of entries to read from the catalogue",
                        type=int, default=-1, action="store")
    parser.add_argument("--mag", 
                        help="largest apparent magnitude to still consider for selecting stars from catalogue",
                        type=float, default=5.5, action="store")
    parser.add_argument("--highlight", dest='highlight', help='what constellations (given by their abbreviation) to hightlight', nargs='+')
    parser.add_argument("--shiftra", 
                        help="shift the RA of all stars by this value (in h)",
                        type=float, default=0, action="store")
    parser.add_argument("--shiftdec", 
                        help="shift the DEC of all stars by this value (in deg)",
                        type=float, default=0, action="store")
    args = parser.parse_args()

    # read in star catalogue data
    starlist = []
    polaris = Star.empty()
    with open('HYG-Database/hygdata_v3.csv', 'r') as csvfile:
        starreader = csv.DictReader(csvfile, delimiter=',')
        for row in starreader:
            s = Star.fromcatalogue(row)
            if s.name == 'Polaris':
                polaris = s
                log.info('Found Polaris in db: {}'.format(s.vec))
            if s.mag < args.mag and s.mag > -26: # filter low-visibility stars and sun
                if args.highlight and row['con'] in args.highlight:
                    log.debug("Found {} with mag {} in constellation '{}' to highlight".format(s.name, s.mag, row['con']))
                    s.highlight = True
                if args.highlight and s.name in args.highlight:
                    log.debug("Found star '{}' to highlight".format(s.name))
                    s.highlight = True
                starlist.append(s)
            if args.nentries > 0 and len(starlist) > args.nentries:
                print ("Reached requested maximum number of events: {}".format(args.nentries))
                break

    log.info("Loaded {} stars into memory".format(len(starlist)))
            
    # plot some of the relevant star data
    plt.figure()
    plt.hist(np.array([s.dec for s in starlist]), 50)
    plt.xlabel("Dec")
    plt.show()
    plt.figure()
    plt.hist(np.array([s.ra for s in starlist]), 50)
    plt.xlabel("RA")
    plt.show()
    plt.figure()
    maglist = np.array([s.mag for s in starlist])
    plt.hist(maglist, 50)
    plt.xlabel("mag")
    plt.show()

    # rotate the stars around axis pointing to polaris eg to align constallations better with face boundaries
    # !!! does not update RA and DEC and therefore might get overwritten !!!
    #M0 = rotation_matrix(polaris.vec, math.pi/4)
    #for s in starlist:
    #    s.vec = np.dot(M0, s.vec)
    if args.shiftra or args.shiftdec:
        log.info('Shifting stars by RA {} and DEC {}'.format(args.shiftra, args.shiftdec))
        for s in starlist:
            s.shift(args.shiftra, args.shiftdec)
        polaris.shift(args.shiftra, args.shiftdec)
    
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
    ax.add_collection3d(Poly3DCollection(verts, facecolor='g'))
    ax.scatter([s.dodecvec[0] for s in pinnedstarlist],[s.dodecvec[1] for s in pinnedstarlist],
               [s.dodecvec[2] for s in pinnedstarlist],
               s=[msize(s.mag) for s in pinnedstarlist],
               color=["red" if s.highlight else "blue" for s in pinnedstarlist])
    plt.show()
    
    # for one face, construct a 2D image of the face on the face's 
    # local xy-plane
    # orientation: local x in global z direction
    for i in range(12):
        # construct unit vectors in the face's plane
        vecface = Dodec.getFaceCtr(i)
        veclocy = np.cross(vecface, np.array([0, 0, 1]))
        veclocy *= 1/math.sqrt(np.dot(veclocy, veclocy))
        veclocx = np.cross(vecface, veclocy)
        veclocx *= 1/math.sqrt(np.dot(veclocx, veclocx))
        # verify vector orientation
        if np.dot(veclocx, [0, 0, 1]) < 0:
            log.debug("Face {} x vector reorientation needed".format(i))
            veclocx *= -1
        if np.dot(np.cross(veclocx, veclocy), vecface) < 0:
            log.debug("Face {} y vector reorientation needed".format(i))
            veclocy *= -1
            
        facestarlist = [s for s in pinnedstarlist if s.faceid == i]
        plt.figure()
        fig, ax = plt.subplots(figsize=(6, 6))
        # draw borders of face
        from matplotlib.patches import Circle, PathPatch
        from matplotlib.path import Path
        vertices = []
        for v in Dodec.getVertices(i):
            norm = 1/np.dot(v, vecface)
            x = np.dot(norm*v-vecface, veclocx)
            y = np.dot(norm*v-vecface, veclocy)
            vertices.append([x,y])
        # close path
        vertices.append(vertices[0])
        # loop over sides and identify bordering faces
        borderfaces = []
        for s in range(5):
            f = _findFacesSharingVertices([Dodec.getVertices(i)[s], 
                                        Dodec.getVertices(i)[(s+1)%5]], Dodec)
            log.debug('Side {} of face {} is shared with face {}'.format(s, i, [_f for _f in f if _f != i]))
            borderfaces.append([_f for _f in f if _f != i][0])
        # add path to draw hexagon outline
        hexagon_path = Path(vertices, closed=True)
        hexagon_patch = PathPatch(hexagon_path, facecolor='none',
                                  edgecolor=(0.0, 0.0, 0.0))
        ax.add_patch(hexagon_patch)
        # plot stars
        plt.scatter([np.dot(s.vec-vecface, veclocx) for s in facestarlist],
                    [np.dot(s.vec-vecface, veclocy) for s in facestarlist],
                    s=[msize(s.mag) for s in facestarlist],
                    color=["red" if s.highlight else "blue" for s in facestarlist])
        # add text indicating bordering faces and this face's index
        for s in range(5):
            textcoord = 1.15*(np.array(vertices[s+1])+0.5*(np.array(vertices[s])-np.array(vertices[s+1])))
            plt.text(textcoord[0], textcoord[1], 'face {}'.format(borderfaces[s]))
        plt.axis((-1, 1, -1, 1))
        plt.show() ## <- shows the plots
    
    ## final commands to show resulting plots (not needed in an interactive session)
    #plt.draw()
    #plt.pause(1) # <-------
    #input("<Hit Enter To Close>")
    #plt.close('all')

if __name__ == "__main__":
	main()
