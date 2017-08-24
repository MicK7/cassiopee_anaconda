"""Initialization of grid solutions.
"""
__version__ = '2.4'
__author__ = "Stephanie Peron, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud"

import initiator

def initConst(array, adim='adim1', MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8):
    """Init the array by a constant field.
    Usage: initConst(array, MInf, alphaZ, alphaY, ReInf)"""
    try:
        import Converter as C
        import KCore.Adim as Adim
    except:
        raise ImportError("initConst: requires Converter module.")
    if adim == 'adim1': cons = Adim.adim1(MInf, alphaZ, alphaY, ReInf)
    else: cons = Adim.adim2(MInf, alphaZ, alphaY, ReInf)

    array = C.initVars(array, 'ro', cons[0])
    array = C.initVars(array, 'rou', cons[1])
    array = C.initVars(array, 'rov', cons[2])
    array = C.initVars(array, 'row', cons[3])
    array = C.initVars(array, 'roE', cons[4])
    return array

def initLamb(array, position=(0.,0.), Gamma=2., MInf=0.5):
    """Init the array defining a grid with a Lamb vortex of
    intensity Gamma and position (x0,y0).
    Usage: initLamb(array, (x0,y0), Gamma, MInf)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(initiator.initLamb(i, position, Gamma, MInf))
        return b
    else:
        return initiator.initLamb(array, position, Gamma, MInf)

def initVisbal(array, position=(0.,0.), Gamma=2., MInf=0.5):
    """Init the array defining a grid with a Visbal vortex of
    intensity Gamma and position (x0,y0).
    Usage: initVisbal(array, (x0,y0), Gamma, MInf)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(initiator.initVisbal(i, position, Gamma, MInf))
        return b
    else:
        return initiator.initVisbal(array, position, Gamma, MInf)

def initYee(array, position=(0.,0.,0.), Gamma=2., MInf=0.5):
    """Init the array defining a grid with a Yee vortex of
    intensity Gamma and position (x0,y0,z0).
    Usage: initYee( array, (x0,y0,z0), Gamma, Minf )"""
    return initiator.initYee( array, position, Gamma, MInf )

def initScully(array, position=(0.,0.), Gamma=2., \
               coreRadius=1., MInf=0.5, model=0):
    """Init the array defining a block field with a Scully vortex
    of intensity Gamma, core radius coreRadius and position (x0,y0).
    Usage:
    initScully(array, (x0,y0), Gamma, coreRadius, MInf, model)"""
    if isinstance(array[0], list):
        b = []
        for i in array:
            b.append(initiator.initScully(i, position, Gamma, coreRadius,
                                          MInf, model))
        return b
    else:
        return initiator.initScully(array, position, Gamma,
                                    coreRadius, MInf, model)

def overlayField(array1, array2, MInf=0.5):
    """Overlay the field of array1 and array2.
    Usage: overlayField(array1, array2, MInf)"""
    if isinstance(array1[0], list):
        b = []; c = 0
        for i in array1:
            b.append(initiator.overlayField(i, array2[c], Gamma, MInf))
            c += 1
        return b
    else:
        return initiator.overlayField(array1, array2, MInf)
