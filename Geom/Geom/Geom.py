"""Geometry definition module.
"""
__version__ = '2.4'
__author__ = "Stephanie Peron, Christophe Benoit, Pascal Raud, Sam Landier"
# 
# Python Interface to define geometries in arrays
#
import geom
import numpy

# - Basic entities -
def point(P):
    """Create a point. Usage: a = point((x,y,z))"""
    a = numpy.zeros((3, 1), numpy.float64)
    a[0,0] = P[0]; a[1,0] = P[1]; a[2,0] = P[2]
    c = numpy.ones((1, 0), numpy.int32)
    return ['x,y,z', a, c, 'NODE']

def naca(epaisseur, N=101):
    """Create a naca profile of N points. Usage: a = naca(eps, N)"""
    return geom.naca(epaisseur, N)    
    
def line(P1, P2, N=100):
    """Create a line of N points. Usage: a = line((x1,y1,z1), (x2,y2,z2), N)"""
    return geom.line(P1, P2, N)

def spline(Pts, order=3, N=100, M=100, density=-1):
    """Create a spline of N points. Usage: a = spline(ctrlsPts, order, N)"""
    return geom.spline(Pts, order, N, order, M, density)

def nurbs(Pts, Weights, order=3, N=100, M=100, density=-1):
    """Create a nurbs of N points. Usage: a = nurbs(ctrlsPts, order, N)"""
    return geom.nurbs(Pts, Weights, order, N, order, M, density)

def polyline(Pts):
    """Create a polyline of N points. Usage: a = polyline([(x1,y1,z1),....,(xn,yn,zn)])"""
    return geom.polyline(Pts)

def circle(C, R, tetas=0., tetae=360., N=100):
    """Create a portion of circle of N points and of center C,
    radius R, between angle tetas and tetae.
    Usage: a = circle((xc,yc,zc), R, tetas, tetae, N)"""
    return geom.circle(C, R, tetas, tetae, N)

def sphere(C, R, N=100):
    """Create a sphere of Nx2N points and of center C and radius R.
    Usage: a = sphere((xc,yc,zc), R, N)"""
    return geom.sphere(C, R, N)

def sphere6(Center, R, N=100):
    """Create a sphere of 6NxN points and of center C and radius R, made of 6 zones.
    Usage: a = sphere6((xc,yc,zc), R, N)"""
    try: import Transform as T; import Generator as G
    except:
        raise ImportError("sphere6: requires Transform and Generator modules.")
    
    s = sphere(Center, R, 2*N)
    b = G.cart((-R/2.+Center[0],-R/2.+Center[1],-R/2.+Center[2]),
               (R/(N-1.),R/(N-1.),R/(N-1.)), (N,N,N))

    b1 = T.subzone(b, (1,1,1), (N,N,1))
    b2 = T.subzone(b, (1,1,1), (1,N,N))
    b3 = T.subzone(b, (1,1,1), (N,1,N))
    b4 = T.subzone(b, (N,1,1), (N,N,N))
    b5 = T.subzone(b, (1,1,N), (N,N,N))
    b6 = T.subzone(b, (1,N,1), (N,N,N))

    c1 = T.projectRay(b1, [s], Center); c1 = T.reorder(c1, (-1,2,3))
    c2 = T.projectRay(b2, [s], Center); c2 = T.reorder(c2, (3,2,1))
    c3 = T.projectRay(b3, [s], Center); c3 = T.reorder(c3, (1,3,2))
    c4 = T.projectRay(b4, [s], Center); c4 = T.reorder(c4, (3,1,2))
    c5 = T.projectRay(b5, [s], Center)
    c6 = T.projectRay(b6, [s], Center); c6 = T.reorder(c6, (-1,3,2))
    return [c1, c2, c3, c4, c5, c6]

def sphereYinYang(C, R, N=100):
    """Create a sphere of center C and radius R made of two overlapping zones.
    Usage: a = sphereYinYang((xc,yc,zc), R, N)"""
    try: import Transform as T
    except: raise ImportError("sphereYinYang: requires Transform module.")
    fringe = 2
    Ni = 4*(N/2)
    a = sphere(C, R, N=Ni)
    a = T.subzone(a, (Ni/4-2,Ni/4-2,1), (3*Ni/4+2,7*Ni/4+2,1))
    b = T.rotate(a, (0,0,0), (0,1,0), 90.)
    b = T.rotate(b, (0,0,0), (0,0,1), 180.)
    return [a, b]

def cone(C, Rb, Rv, H, N=100):
    """Create a cone of NxNh points and of center C, basis radius Rb, vertex radius Rv and height H.
    Usage: a = cone((xc,yc,zc), Rb, Rv, H, N)"""
    return geom.cone(C, Rb, Rv, H, N)

def torus(C, R, r, alphas=0., alphae=360.,
          betas=0., betae=360., NR=100, Nr=100):
    """Create NRxNr points lying on a torus of center C and radii R (main)
    and r (tube) between the angles alphas and alphae (XY-plane) and between
    betas and betae (RZ-plane).
    Usage: a = torus((xc,yc,zc), R, r, alphas, alphae, betas, betae, NR, Nr)"""
    return geom.torus(C, R, r, alphas, alphae, betas, betae, NR, Nr)

def triangle(P1, P2, P3):
    """Create a single triangle with points P1, P2, P3.
    Usage: a = triangle((x1,y,1,z1), (x2,y2,z2), (x3,y3,z3))"""
    return geom.triangle(P1, P2, P3)

def quadrangle(P1, P2, P3, P4):
    """Create a single quadrangle with points P1, P2, P3, P4.
    Usage: a = quadrangle((x1,y,1,z1), (x2,y2,z2), (x3,y3,z3), (x4,y4,z4))"""
    return geom.quadrangle(P1, P2, P3, P4)

def bezier(controlPts, N=100, M=100, density=-1):
    """Create a a Bezier curve defined by an array of control points controlPts.
    Usage: a = bezier(controlPts, N, M)"""
    return geom.bezier(controlPts, N, M, density)

def curve(f, N=100):
    """Create a curve from a user defined parametric function or a formula.
    Usage: a = curve(f, N)"""
    if isinstance(f, str): return curve_(f, N)
    else: return curve__(f, N)

# Courbe parametree a partir d'une formule
def curve_(f, N):
    import Converter; import Generator
    a = Generator.cart( (0,0,0), (1./(N-1),1,1), (N,1,1))
    a[0] = 't,y,z'
    a = Converter.initVars(a, f)
    a = Converter.extractVars(a, ['x','y','z'])
    return a

# Courbe parametree a partir d'une fonction
def curve__(f, N):
    a = numpy.zeros((3, N), dtype=numpy.float64)
    r = f(0)
    if len(r) != 3:
        print "Warning: curve: parametric function must return a (x,y,z) tuple."
        return ['x,y,z', a, N, 1, 1]
    for i in range(N):
        t = 1.*i/(N-1)
        r = f(t)
        a[0,i] = r[0]; a[1,i] = r[1]; a[2,i] = r[2]
    return ['x,y,z', a, N, 1, 1]

def surface(f, N=100):
    """Create a surface from a user defined parametric function or a formula.
    Usage: a = surface(f, N)"""
    if isinstance(f, str): return surface_(f, N)
    else: return surface__(f, N)

# Surface parametree a partir d'une formule
def surface_(f, N):
    import Converter; import Generator
    a = Generator.cart( (0,0,0), (1./(N-1),1./(N-1),1), (N,N,1))
    a[0] = 't,u,z'
    a = Converter.initVars(a, f)
    a = Converter.extractVars(a, ['x','y','z'])
    return a
    
# Surface parametree a partir d'une fonction
def surface__(f, N):
    a = numpy.zeros((3, N*N), dtype=numpy.float64)
    r = f(0,0)
    if len(r) != 3:
        print "Warning: surface: parametric function must return a (x,y,z) tuple."
        return ['x,y,z', a, N, N, 1]
    for j in range(N):
        u = 1.*j/(N-1)
        for i in range(N):
            ind = i + j*N
            t = 1.*i/(N-1)
            r = f(t,u)
            a[0,ind] = r[0]
            a[1,ind] = r[1]
            a[2,ind] = r[2]
    return ['x,y,z', a, N, N, 1]

# - informations -
def getLength(array):
    """Return the length of 1D array(s) defining a mesh.
    Usage: l = getLength(array(s))"""
    if isinstance(array[0], list): 
        l = 0.
        for i in array: l += geom.getLength(i)
        return l
    else: return geom.getLength(array)

def getDistantIndex(array, ind, l):
    """Return the index of 1D array defining a mesh located at a
    distance l of ind.
    Usage: ind = getDistantIndex(array, ind, l)"""
    return geom.getDistantIndex(array, ind, l)

def getNearestPointIndex(array, pointList):
    """Return the nearest index of points in array.
    Usage: getNearestPointIndex(array, pointList)"""
    if isinstance(pointList, tuple): pL = [pointList]
    else: pL = pointList

    if isinstance(array[0], list):
        # keep nearest
        npts = len(pL)
        res0 = [(0,1.e6) for i in xrange(npts)]
        for i in array:
            res = geom.getNearestPointIndex(i, pL)
            for j in xrange(npts):
                if (res0[j][1] > res[j][1]):
                    res0[j] = (res[j][0], res[j][1])
        if isinstance(pointList, tuple): return res0[0]
        else: return res0
    else:
        res = geom.getNearestPointIndex(array, pL)
        if isinstance(pointList, tuple): return res[0]
        else: return res

def getCurvatureRadius(array):
    """Return the curvature radius for each point.
    Usage: getCurvatureRadius(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(geom.getCurvatureRadius( i ))
        return b
    else:
        return geom.getCurvatureRadius(array)
    
def getCurvatureAngle(array):
    """Return the curvature angle for each point...
    Usage: getCurvatureAngle(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(geom.getCurvatureAngle(i))
        return b
    else:
        return geom.getCurvatureAngle(array)
    
def getCurvatureHeight(array):
    """Return the curvature height for each node in a 2D or 1D array...
    Usage: getCurvatureHeight(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(geom.getCurvatureHeight(i))
        return b
    else:
        return geom.getCurvatureHeight(array)

def getSharpestAngle(array):
    """Return the sharpest angle for each point of a surface based on the sharpest angle
    between adjacent element to which the point belongs to.
    Usage: getSharpestAngle(a)"""
    if isinstance(array[0], list): 
        out = []
        for i in array: out.append(geom.getSharpestAngle(i))
        return out
    else: return geom.getSharpestAngle(array)
    
def getCurvilinearAbscissa(array):
    """Return the curvilinear abscissa for each point...
    Usage: getCurvilinearAbscissa(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(geom.getCurvilinearAbscissa(i))
        return b
    else:
        return geom.getCurvilinearAbscissa(array)
        
def getDistribution(array):
    """Return the curvilinear abscissa for each point as coordinates
    Usage: getDistribution(array)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            if i[-1]=='BAR': raise TypeError("getDistribution: only for structured array.")
            c = line((0,0,0),(1,0,0),i[2])
            c[1][0] = geom.getCurvilinearAbscissa(i)[1]
            b.append(c)
        return b
    else:
        if array[-1]=='BAR': raise TypeError("getDistribution: only for structured arrays.")
        c = line((0,0,0),(1,0,0),array[2])
        c[1][0] = geom.getCurvilinearAbscissa(array)[1]    
        return c       
        
def getTangent(array):
    """Makes the tangent of a 1D curve, as a new array. The input argument 
    shall be an array. Each node of the output represents the unitary tangent 
    vector, pointing towards the tangent direction of the input 1D curve.
    Usage: getTangent(array)"""
    np = numpy
    if not isinstance(array[0], list): Arrays = [array]
    else: Arrays = array
    b = []
    for a in Arrays:
        t = ['x,y,z',0,a[2],1,1]
        # Central difference
        n = a[1]
        OrientationAbsolute = 0.5*(np.diff(n[:,:-1],axis=1)+np.diff(n[:,1:],axis=1))
        # Not centered for the bounds
        OrientationAbsolute = np.hstack(((n[:,1]-n[:,0])[np.newaxis].T,
                                         OrientationAbsolute,
                                         (n[:,-1]-n[:,-2])[np.newaxis].T))
        #Norm = np.linalg.norm(OrientationAbsolute, axis=0)
        Norm = np.sqrt(np.sum(OrientationAbsolute*OrientationAbsolute, axis=0))
        OrientationRelative = OrientationAbsolute/Norm
        t[1] = OrientationRelative
        b.append(t)
    if len(b)==1: return b[0]
    else: return b

def lineGenerate(array, arrayLine):
    """Generate a surface mesh by using 1D array (defining a mesh)
    and following the curve defined in arrayLine.
    Usage: lineGenerate(array, arrayLine)"""
    if isinstance(arrayLine[0], list): # set of driving curves
        if isinstance(array[0], list):
            b = []
            for i in array:
                b.append(lineGenerate2__(i, arrayLine))
            return b
        else:
            return lineGenerate2__(array, arrayLine)
    else: # one driving curve
        if isinstance(array[0], list):
            b = []
            for i in array:
                b.append(geom.lineGenerate(i, arrayLine))
            return b
        else:
            return geom.lineGenerate(array, arrayLine)

def lineGenerate2__(array, drivingCurves):
    import Converter; import Generator
    # Copie la distribution de 0 sur les autres courbes
    d = []
    ref = drivingCurves[0]; d += [ref]
    l = getLength(ref)
    distrib = getCurvilinearAbscissa(ref)
    distrib[0] = 'x'; distrib = Converter.addVars(distrib, ['y','z'])
    for i in drivingCurves[1:]:
        d += [Generator.map(i, distrib)]
    return geom.lineGenerate2(array, d)

def addSeparationLine(array, array2):
    """Add a separation line defined in array2 to a mesh defined in array.
    Usage: addSeparationLine(array, array2)"""
    return geom.addSeparationLine(array, array2)

def axisym(array, center, axis, angle=360., Ntheta=180, rmod=None):
    """Create an axisymetriccal mesh given an azimuthal surface mesh.
    Usage: axisym(array, (xo,yo,zo), (nx,ny,nz), teta, Nteta, rmod)"""
    try: 
        import Converter
        if rmod is not None: rmod = Converter.convertBAR2Struct(rmod)
    except: pass 
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(geom.axisym(i, center, axis, angle, Ntheta, rmod))
        return b
    else:
        return geom.axisym(array, center, axis, angle, Ntheta, rmod)

def volumeFromCrossSections(contours):
    """Generate a 3D volume from cross sections contours in (x,y) planes.
    Usage: volumeFromCrossSections(contours)"""
    try:
        import KCore
        import Converter as C
        import Transform as T
        import Generator as G
    except:
        raise ImportError("volumeFromCrossSections: require Converter, Transform and Generator.")

    c = {}
    # Dictionnaire des contours suivant z
    for i in contours:
        posz = KCore.isNamePresent(i, "z")
        if posz == -1:
            posz = KCore.isNamePresent(i, "Z")
            if posz == -1:
                posz = KCore.isNamePresent(i, "CoordinateZ")
                if posz == 1:
                    raise TypeError("volumeFromCrossSections: Z coordinates not found in an array.")
        z = C.getValue(i, 0)[posz]
        if z in c:
            d = C.convertArray2Tetra(i)
            b = c[z]; f = T.join(d, b); c[z] = f
        else:
            d = C.convertArray2Tetra(i)
            c[z] = C.convertArray2Tetra(d)

    sort = sorted(c.items())
    
    # Delaunay constrained
    DT = []; CT = []
    for i in sort:
        d = G.close(i[1], 1.e-6); CT.append(d)
        #m = G.constrainedDelaunay(d, 1.e-10, 1)
        m = G.T3mesher2D(d); m = T.translate(m, (0,0,i[0]))
        DT.append(m)
    if len(sort) < 2:
        raise ValueError("volumeFromCrossSections: require at least two cross sections.")

    l = 0
    for i in sort:
        if l == 0:
            vol = geom.volumeFromCrossSections(DT[l], DT[l+1], CT[l], CT[l+1])
        elif l < len(sort)-1:
            vol2 = geom.volumeFromCrossSections(DT[l], DT[l+1], CT[l], CT[l+1])
            vol = T.join(vol, vol2)
        l += 1
    return vol

# - text functions -
def text1D(string, font='text1', smooth=0, offset=0.5):
    """Create a 1D text. offset is the space between letters.
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text1D(string, font, smooth, offset)"""
    if font == 'text1': import text1 as Text
    elif font == 'vera': import vera as Text
    elif font == 'chancery': import chancery as Text
    elif font == 'courier': import courier as Text
    elif font == 'nimbus': import nimbus as Text
    else: import text1 as Text
    try: import Transform
    except: raise ImportError("text1D: requires Transform.")
    retour = []
    offx = 0.; offy = 0.; s = 6
    for i in string:
        if (i == 'A'): a, s = Text.A()
        elif (i == 'a'): a, s = Text.a()
        elif (i == 'B'): a, s = Text.B()
        elif (i == 'b'): a, s = Text.b()
        elif (i == 'C'): a, s = Text.C()
        elif (i == 'c'): a, s = Text.c()
        elif (i == 'D'): a, s = Text.D()
        elif (i == 'd'): a, s = Text.d()
        elif (i == 'E'): a, s = Text.E()
        elif (i == 'e'): a, s = Text.e()
        elif (i == 'F'): a, s = Text.F()
        elif (i == 'f'): a, s = Text.f()
        elif (i == 'G'): a, s = Text.G()
        elif (i == 'g'): a, s = Text.g()
        elif (i == 'H'): a, s = Text.H()
        elif (i == 'h'): a, s = Text.h()
        elif (i == 'I'): a, s = Text.I()
        elif (i == 'i'): a, s = Text.i()
        elif (i == 'J'): a, s = Text.J()
        elif (i == 'j'): a, s = Text.j()
        elif (i == 'K'): a, s = Text.K()
        elif (i == 'k'): a, s = Text.k()
        elif (i == 'L'): a, s = Text.L()
        elif (i == 'l'): a, s = Text.l()
        elif (i == 'M'): a, s = Text.M()
        elif (i == 'm'): a, s = Text.m()
        elif (i == 'N'): a, s = Text.N()
        elif (i == 'n'): a, s = Text.n()
        elif (i == 'O'): a, s = Text.O()
        elif (i == 'o'): a, s = Text.o()
        elif (i == 'P'): a, s = Text.P()
        elif (i == 'p'): a, s = Text.p()
        elif (i == 'Q'): a, s = Text.Q()
        elif (i == 'q'): a, s = Text.q()
        elif (i == 'R'): a, s = Text.R()
        elif (i == 'r'): a, s = Text.r()
        elif (i == 'S'): a, s = Text.S()
        elif (i == 's'): a, s = Text.s()
        elif (i == 'T'): a, s = Text.T()
        elif (i == 't'): a, s = Text.t()
        elif (i == 'U'): a, s = Text.U()
        elif (i == 'u'): a, s = Text.u()
        elif (i == 'V'): a, s = Text.V()
        elif (i == 'v'): a, s = Text.v()
        elif (i == 'W'): a, s = Text.W()
        elif (i == 'w'): a, s = Text.w()
        elif (i == 'X'): a, s = Text.X()
        elif (i == 'x'): a, s = Text.x()
        elif (i == 'Y'): a, s = Text.Y()
        elif (i == 'y'): a, s = Text.y()
        elif (i == 'Z'): a, s = Text.Z()
        elif (i == 'z'): a, s = Text.z()
        elif (i == '0'): a, s = Text.C0()
        elif (i == '1'): a, s = Text.C1()
        elif (i == '2'): a, s = Text.C2()
        elif (i == '3'): a, s = Text.C3()
        elif (i == '4'): a, s = Text.C4()
        elif (i == '5'): a, s = Text.C5()
        elif (i == '6'): a, s = Text.C6()
        elif (i == '7'): a, s = Text.C7()
        elif (i == '8'): a, s = Text.C8()
        elif (i == '9'): a, s = Text.C9()
        elif (i == '.'): a, s = Text.POINT()
        elif (i == ','): a, s = Text.COMMA()
        elif (i == ';'): a, s = Text.POINTCOMMA()
        elif (i == ':'): a, s = Text.TWOPOINTS()
        elif (i == '!'): a, s = Text.EXCLAMATION()
        elif (i == '+'): a, s = Text.PLUS()
        elif (i == '-'): a, s = Text.MINUS()
        elif (i == '='): a, s = Text.EQUAL()
        elif (i == '('): a, s = Text.LEFTBRACE()
        elif (i == ')'): a, s = Text.RIGHTBRACE()
#        elif (i == u'\xe9'):
#            a, s = Text.EACUTE()
        elif (i == '\n'):
            offy = offy - 8 - offset
            offx = -6 - offset
            a = []
        else:
            a = []
        if a != []:
            a = Transform.translate(a, (offx,offy,0))
            retour += a
        offx += s + offset

    if smooth != 0:
        try: import Generator
        except:
            raise ImportError("text1D: requires Generator for smooth option.")
        if smooth == 1:
            nmap = 40; hdensify = 8./100
        elif smooth == 2:
            nmap = 40; hdensify = 8./10.
        elif smooth == 3:
            nmap = 40; hdensify = 8./5
        else:
            nmap = 40; hdensify = 8./2
        d = Generator.cart((0,0,0), (1./nmap,1,1), (nmap+1,1,1))
        c = 0
        for i in retour:
            b = Generator.densify(i, hdensify) 
            retour[c] = Generator.map(b, d); c += 1
            
    return retour

def text2D(string, font='text1', smooth=0, offset=0.5):
    """Create a 2D text. offset is the space between letters.
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text2D(string, font, smooth, offset)"""
    try:
        import Generator; import Transform; import Converter
    except:
        raise ImportError("text2D: requires Generator, Transform, Converter.")
    a = text1D(string, font, smooth, offset)
    a = Converter.convertArray2Tetra(a)
    b = Transform.join(a)
    b = Generator.constrainedDelaunay(b)
    b = Generator.selectInsideElts(b, a)
    return b

def text3D(string, font='text1', smooth=0, offset=0.5):
    """Create a 3D text. offset is the space between letters.
    font indicates used font, smooth indicates the smoothing intensity.
    Usage: text3D(string, font, smooth, offset)"""
    try:
        import Generator
        import Transform
        import Converter
    except:
        raise ImportError("text3D: requires Generator, Transform, Converter.")
    a = text1D(string, font, smooth, offset)
    l = line((0,0,0),(0,0,8),2)
    a = lineGenerate(a, l)
    a = Converter.convertArray2Tetra(a)
    b = Transform.join(a)
    
    a = text2D(string, font, smooth, offset)
    a = Transform.translate(a, (0,0,-0.0001))
    c = Transform.translate(a, (0,0,8.0002))
    a = Transform.join([a, b, c])
    return a

#======================================================================
# connect 1D curves
# IN: sharpness=0: par des lignes, =1 par des splines
# IN: N: nbre de pts dans les raccords
# IN: lengthFactor: enleve les raccords trop longs
#======================================================================
def connect1D(curves, sharpness=0, N=10, lengthFactor=1.):
    """Connect 1D curves in a single curve.
    Usage: a = connect1D(A, sharpness, N, lengthFactor)"""
    import KCore.Vector as Vector
    import Transform as T
    import Converter as C
    import Generator as G    
    #curves = T.splitTBranch(curves)
    curves = C.convertBAR2Struct(curves)
    ncurves = len(curves)

    Pts = []; PtsM= []
    for i in curves:
        ni = i[2]
        e1 = C.getValue(i, 0)
        e2 = C.getValue(i, ni-1)
        e1M = C.getValue(i, 1)
        e2M = C.getValue(i, ni-2)
        Pts.append([e1, e2])
        PtsM.append([e1M, e2M])

    added = []
    for c in xrange(ncurves):
        lcurve = getLength(curves[c]) * lengthFactor
        P1 = Pts[c][0] 
        minDist, P2, d, ext = findNearest__(P1, Pts, c)
        n1 = Vector.sub(P1, PtsM[c][0])
        n1 = Vector.normalize(n1)
        n2 = Vector.sub(P2, PtsM[d][ext]) 
        PI = intersectionPoint__(P1,n1,P2,n2)
        if sharpness == 0: # sharp
            la = line(P1,PI, N=N)
            lb = line(PI,P2, N=N)
            if getLength(la) < lcurve and getLength(lb) < lcurve:
                added += [la,lb]
        elif sharpness == 1: # spline
            controlPts = polyline([P1,PI,P2])
            sp = spline(controlPts, N=N)
            if getLength(sp) < lcurve: added += [sp]

        P1 = Pts[c][1]
        minDist, P2, d, ext = findNearest__(P1, Pts, c)
        n1 = Vector.sub(P1, PtsM[c][0])
        n1 = Vector.normalize(n1)
        n2 = Vector.sub(P2, PtsM[d][ext])
        n2 = Vector.normalize(n2)
        PI = intersectionPoint__(P1,n1,P2,n2)
        if sharpness == 0: # sharp
            la = line(P1,PI, N=N)
            lb = line(PI,P2, N=N)
            if getLength(la) < lcurve and getLength(lb) < lcurve:
                added += [la,lb]
        elif sharpness == 1: # spline
            controlPts = polyline([P1,PI,P2])
            sp = spline(controlPts, N=N)
            if getLength(sp) < lcurve: added += [sp]
    out = C.convertArray2Hexa(curves)
    if added != []:
        added = C.convertArray2Hexa(added)
        out = T.join(out+added)
    out = G.close(out) 
    return out

# Pt d'intersection par minimal distance
def intersectionPoint__(P1,n1,P2,n2):
    import KCore.Vector as Vector
    s = Vector.dot(n1, n2)
    s2 = 1.-s*s
    dP1P2 = Vector.sub(P1, P2)
    sn1 = Vector.mul(s, n1)
    p2 = Vector.sub(n2, sn1)
    q = Vector.dot(dP1P2, p2)
    tp = q / s2
    t = Vector.dot(Vector.sub(P2,P1),n1)+s*tp
    PI1 = Vector.add(P2, Vector.mul(tp,n2))
    PI2 = Vector.add(P1, Vector.mul(t, n1))
    PI = Vector.add(PI1,PI2)
    PI = Vector.mul(0.5,PI)
    return PI

# trouve le pt le plus proche de Pt dans Pts mais different de c
def findNearest__(Pt, Pts, c):
    import KCore.Vector as Vector
    minDist = 1.e6; nearest = None; dmin = -1; ext=0;
    for d in xrange(len(Pts)):
        if d <= c: # possible sur lui meme !!
            e2a = Pts[d][0]; e2b = Pts[d][1]    
            d1 = Vector.squareDist(Pt, e2a)
            d2 = Vector.squareDist(Pt, e2b)
            if d1 < minDist and d1 > 1.e-12: 
                minDist = d1; nearest = e2a; ext=0; dmin = d
            if d2 < minDist and d2 > 1.e-12:
                minDist = d2; nearest = e2b; ext=1; dmin = d
    return minDist, nearest, dmin, ext

#=======================================================================
# Ferme tous les trous dans a
#=======================================================================
def closeHoles(a):
    ext = P.exteriorFaces(a)
    ext = C.convertArray2Tetra(ext)
    ext = T.join(ext)
    ext = G.close(ext)
    ext = T.splitConnexity(ext)
    out = []
    for i in ext:
            p = G.fittingPlaster(i, bumpFactor=0.)
            b = G.gapfixer(i, p)
            out.append(b)
    return out 
