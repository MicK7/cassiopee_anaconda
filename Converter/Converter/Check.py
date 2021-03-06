# - check -
# Un module de verification de la coherence des arbres pythons

import Internal
import numpy
import PyTree as C

# last is 76
CGNSTypes = {
    'CGNSTree_t':72,
    'CGNSLibraryVersion_t':1,
    'CGNSBase_t':0,
    'BaseIterativeData_t':20,
    'Zone_t':48,
    'Elements_t':52,
    'GridCoordinates_t':54,
    'FlowSolution_t':53,
    'ZoneGridConnectivity_t':64,
    'GridConnectivityProperty_t':65,
    'GridConnectivityType_t':66,
    'GridConnectivity1to1_t':67,
    'GridConnectivity_t':75,
    'Periodic_t':76,
    'OversetHoles_t':68,
    'ZoneIterativeData_t':69,
    'ZoneSubRegion_t':70,
    'ZoneType_t':71,

    'DataArray_t':3,
    'DataConversion_t':4,
    'DataClass_t':5,
    'Descriptor_t':6,
    'DimensionalExponents_t':7,
    'AdditionalExponents_t':8,
    'DimensionalUnits_t':9,
    'AdditionalUnits_t':10,
    'UserDefinedData_t':13,
    'DiscreteData_t':51,
    'Rind_t':50,
    'GridLocation_t':16,
    'Ordinal_t':17,
    'IndexArray_t':18,
    'IndexRange_t':19,
    '"int[IndexDimension]"':73,
    '"int"':74,

    'Axisymmetry_t':2,
    'AxisymmetryAxisVector_t':11,
    'AxisymmetryReferencePoint_t':12,

    'FamilyName_t':15,
    'AdditionalFamilyName_t':14,
    'Family_t':21,
    'FamilyBC_t':22,
    'FamilyBCDataSet_t':23,
        
    'GeometryReference_t':27,
    'GeometryEntity_t':28,
    'GeometryFile_t':29,
    'GeometryFormat_t':30,

    'ReferenceState_t':25,
    'RotatingCoordinates_t':31,
    'FlowEquationSet_t':32,
    'ChemicalKineticsModel_t':33,
    'EMConductivityModel_t':34,
    'EMElectricFieldModel_t':35,
    'GasModel_t':36,
    'GoverningEquations_t':37,
    'ThermalConductivityModel_t':38,
    'ThermalRelaxationModel_t':39,
    'TurbulenceClosure_t':40,
    'TurbulenceModel_t':41,
    'ViscosityModel_t':42,
    'ConvergenceHistory_t':43,
    'Gravity_t':44,
    'IntegralData_t':45,
    'ReferenceState_t':46,
    'SimulationType_t':47,
  
    'ArbitraryGridMotion_t':49,
    'RigidGridMotion_t':55,
    
    'ZoneBC_t':56,
    'BC_t':57,
    'BCData_t':24,
    'BCDataSet_t':58,
    'BCProperty_t':59,
    'Area_t':60,
    'AreaType_t':61,
    'WallFunction_t':62,
    'WallFunctionType_t':63
}

#==============================================================================
# IN: t: pyTree to be checked
# IN: level: check level 0 (version node), 1 (node conformity),
# 2 (unique base name), 3 (unique zone name), 4 (unique BC name),
# 5 (BC ranges), 6 (BCMatch/NearMatch), 7 (FamilyZone et FamilyBCs),
# 8 (invalid CGNS Types), 9 (invalid connectivity)
# if level=-n, perform check from 0 to n
#==============================================================================
def checkPyTree(t, level=-20):
    errors = []
    if (level <= 0 or level == 0):
        # check version node
        errors += checkVersionNode(t)
    if (level <= -1 or level == 1):
        # check nodes conformity
        errors += checkNodes(t)
    if (level <= -2 or level == 2):
        # check unique base names
        errors += checkUniqueNames(t, 'CGNSBase_t')
    if (level <= -3 or level == 3):
        # check unique zone names
        errors += checkUniqueNames(t, 'Zone_t')
    if (level <= -4 or level == 4):
        # check unique BC names
        errors += checkUniqueNames(t, 'BC_t')
        # check unique BCMatch names
        errors += checkUniqueNames(t, 'GridConnectivity1to1_t')
        # check unique BCNearMatch/BCOverlap names
        errors += checkUniqueNames(t, 'GridConnectivity_t')
    if (level <= -5 or level == 5):
        # check BC range
        errors += checkBCRanges(t, 'BC_t')
        # check BCMatch range
        errors += checkBCRanges(t, 'GridConnectivity1to1_t')
        # check BCMatch opposite range
        errors += checkDonorRanges(t, 'GridConnectivity1to1_t')
        # check BCNearMatch/Overlap range
        errors += checkBCRanges(t, 'GridConnectivity_t')
        # check BCNearMatch opposite range
        errors += checkDonorRanges(t, 'GridConnectivity_t')
    if (level <= -6 or level == 6):
        # check BCMatch opposite range
        errors += checkOppositRanges(t, 'GridConnectivity1to1_t')
        # check BCNearMatch/Overlap opposite range
        errors += checkOppositRanges(t, 'GridConnectivity_t')
    if (level <= -7 or level == 7):
        # check zone family
        errors += checkZoneFamily(t)
        # check BC family
        errors += checkBCFamily(t)
    if (level <= -8 or level == 8):
        # check CGNS type for each node
        errors += checkCGNSType(t)
    if (level <= -9 or level == 9):
        # check element nodes (connectivity)
        errors += checkElementNodes(t)
    if (level <= -10 or level == 10):
        # check valid CGNS var name
        errors += checkCGNSVarNames(t)
    return errors

#==============================================================================
# Correct pyTree
#==============================================================================
def correctPyTree(t, level=-20):
    tp = Internal.copyRef(t)
    _correctPyTree(tp, level)
    return tp

def _correctPyTree(t, level=-20):
    # Corrige le noeud version
    if (level <= 0 or level == 0):
        _correctVersionNode(t)
    # Supprime les noeuds non conformes
    if (level <= -1 or level == 1):
       _correctNodes(t)
    # Renomme les bases
    if (level <= -2 or level == 2):
        _correctNames(t, 'CGNSBase_t')
    # Renomme les zones
    if (level <= -3 or level == 3):
        _correctNames(t, 'Zone_t')
    # Renomme les BCs
    if (level <= -4 or level == 4):
        _correctNames(t, 'BC_t')
        _correctNames(t, 'GridConnectivity1to1_t')
        _correctNames(t, 'GridConnectivity_t')
    # Supprime les BCs avec des ranges invalides
    if (level <= -5 or level == 5):
        _correctBCRanges(t, 'BC_t')
        _correctBCRanges(t, 'GridConnectivity1to1_t')
        _correctDonorRanges(t, 'GridConnectivity1to1_t')
        _correctBCRanges(t, 'GridConnectivity_t')
        _correctDonorRanges(t, 'GridConnectivity_t')
    # Supprime les BCs avec des fenetres opposees invalides
    if (level <= -6 or level == 6):
        _correctOppositRanges(t, 'GridConnectivity1to1_t')
        _correctOppositRanges(t, 'GridConnectivity_t')
    # Corrige les family oubliees
    if (level <= -7 or level == 7):
        _correctZoneFamily(t)
        _correctBCFamily(t)
    # Corrige les noeud pas de bon type CGNS
    if (level <= -8 or level == 8):
        _correctCGNSType(t)
    # Corrige les noeuds elements (connectivity)
    if (level <= -9 or level == 9):
        _correctElementNodes(t)
    # Corrige les noms de variables non CGNS
    if (level <= -10 or level == 10):
        _correctCGNSVarNames(t)
    C.registerAllNames(t)
    return None

#==============================================================================
# Check version node
# Doit etre en premier dans l'arbre et non duplique.
#==============================================================================
def checkVersionNode(t):
    errors = []
    if (len(t) == 4 and t[3] == 'CGNSTree_t'):
        version = Internal.getNodesFromType1(t, 'CGNSLibraryVersion_t')
        if version == []: errors += [None, 'Missing CGNS version node.']
        elif (t[2][0] != version[0]): errors += [version[0], 'CGNS version node must be first in tree sons.']
        if len(version) > 1:
            for v in version[1:]: errors += [v, 'Only one CGNS version node is allowed.']
    return errors

#==============================================================================
# Correct version node
#==============================================================================
def _correctVersionNode(t):
    errors = checkVersionNode(t)
    le = len(errors)/2
    added = 0
    for e in xrange(le):
        node = errors[2*e]
        if node is None: t[2].insert(0, Internal.createCGNSVersionNode())
        else:
            c = 0
            for n in t[2]:
                if (n[3] == 'CGNSLibraryVersion_t' and id(n) == id(node)):
                    del t[2][c]; break
                c += 1
            if added == 0: t[2].insert(0, node); added = 1 
    return None

#==============================================================================
# Check nodes renvoie la liste des noeuds non conformes
# Check nodes modifie aussi les noeuds ayant comme valeur des chaines
# de caracteres avec des blancs au debut ou a la fin
#==============================================================================
def checkNodes(node):
    errors = []
    isStd = Internal.isStdNode(node)
    if (isStd >= 0):
        for c in node[isStd:]: checkNode__(c, errors)
    else: checkNode__(node, errors)
    return errors

#==============================================================================
def checkNode__(node, errors):
    sons = []
    # node doit etre une liste
    if isinstance(node, list):
        # node doit avoir pour longueur 4 [nom, valeur, [fils], type]
        if len(node) == 4:
            # si node[1] est une string -> strip
            if isinstance(node[1], str): node[1] = node[1].strip()
            if (isinstance(node[1], numpy.ndarray) and (node[1].dtype.kind == 'S' or node[1].dtype.kind == 'a')):
                val = node[1].tostring(); val = val.strip()
                node[1] = numpy.fromstring(val, 'c')
            
            # node[0] (nom) est un string ou None
            if isinstance(node[0], str) == False and node[0] is None:
                errors += [node, "Node[0] of node %s must be a string designing node name."%node[0]]
                
            # node[2] (liste des fils) doit etre une liste
            if not isinstance(node[2], list):
                errors += [node, "Node[2] of node %s must be a list of sons."%node[0]]
            else: sons = node[2]
            
            # node[3] doit etre un string se terminant par _t ...
            if not isinstance(node[3], str):
                errors += [node, "Node[3] of node %s must be a string designing the node type."%node[0]]
            #if node[3].find('_t') != len(node[3])-2:
            #    errors += [node, "Node[3] of node %s must be a string designing the node type."%node[0]]
                   
        else: errors += [node, "Node %s has a length != 4."%node[0]]
    else: errors += [node, "Node is not a list."]

    for n in sons: checkNode__(n, errors)

#==============================================================================
# Delete les noeuds non conformes
#==============================================================================
def _correctNodes(t):
    errors = checkNodes(t)
    le = len(errors)/2
    for e in xrange(le):
        node = errors[2*e]
        (p, c) = Internal.getParentOfNode(t, node)
        if p is not None: del p[2][c]
    return None

#==============================================================================
# Verifie que les noms des noeuds de type donne sont bien uniques
#==============================================================================
def checkUniqueNames(t, type):
    nameServer = {}
    errors = []
    # Register nodes of type
    nodes = Internal.getNodesFromType(t, type)
    for n in nodes:
        name = n[0]
        if not nameServer.has_key(name): nameServer[name] = 0
        else:
            if type == 'CGNSBase_t':
                errors += [n, "Base name %s is already used."%n[0]]
            elif type == 'Zone_t':
                (p,c) = Internal.getParentOfNode(t, n)
                errors += [n, "Zone name %s (base %s) is already used."%(n[0],p[0])]
            elif type == 'BC_t':
                (p,c) = Internal.getParentOfNode(t, n)
                (p,c) = Internal.getParentOfNode(t, p)
                errors += [n, "BC name %s (zone %s) is already used."%(n[0],p[0])]
            elif type == 'GridConnectivity_1to1_t':
                errors += [n, "BCMatch name %s is already used."%n[0]]
            elif type == 'GridConnectivity_t':
                errors += [n, "GridConnectivity name %s is already used."%n[0]]
            else: errors += [n, "Node name %s is already used."%n[0]]
    return errors

#==============================================================================
# Change node name if already defined
#==============================================================================
def _correctNames(t, type):
    nameServer = {}
    nodes = Internal.getNodesFromType(t, type)
    zoneDonors = []
    for n in nodes:
        name = n[0]
        if not nameServer.has_key(name): nameServer[name] = 0
        else: # deja existant
            c = nameServer[name]; ret = 1
            while (ret == 1):
                name2 = '%s.%d'%(name,c)
                if not nameServer.has_key(name2): ret = 0
                else: ret = 1
                c += 1
            nameServer[name2] = 0
            nameServer[name] = c
            if n[3] == 'Zone_t': zoneDonors.append((n[0], name2))
            n[0] = name2
            
    # Modifie les zoneDonors
    _correctDonors(t, 'GridConnectivity1to1_t', zoneDonors)
    _correctDonors(t, 'GridConnectivity_t', zoneDonors)
    
    # Modifie les attachs
    return None

#==============================================================================
def _correctDonors(t, type, zoneDonors):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType2(z, type)
        for n in nodes:
            if isinstance(n[1], numpy.ndarray):
                zdonorname = n[1].tostring()
            else: zdonorname = n[1]
            for zd in zoneDonors:
                zd1 = numpy.fromstring(zd[1], 'c')
                if zd[0] == zdonorname: n[1] = zd1
    return None

#==============================================================================
# Check BC ranges
# IN: t: arbre a verifier
# IN: type: type du noeud de BC a verifier (BC_t,...)
# Verifie que le range de la BC est contenu dans la grille
# Verifie que les faces de la BC sont contenues dans la grille
# Corrige la shape des ranges si celle-ci est en C
#==============================================================================
def checkBCRanges(t, type):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        nodes = Internal.getNodesFromType2(z, type)
        for n in nodes:
            prange = Internal.getNodesFromName1(n, 'PointRange')
            for r in prange:
                if (r[1].shape == (2,3)):
                    r[1] = numpy.reshape(r[1], (3,2), order='Fortran')
                if (r[1].shape == (3,2)):
                    win = Internal.range2Window(r[1])
                else: win = [0,0,0,0,0,0] # pas de check en non structure
                # Check structure uniquement pour l'instant
                error = 0
                if (win[0] < 0 or win[0] > dim[1]): error = 1
                if (win[1] < 0 or win[1] > dim[1]): error = 1
                if (win[2] < 0 or win[2] > dim[2]): error = 1
                if (win[3] < 0 or win[3] > dim[2]): error = 1
                if (win[4] < 0 or win[4] > dim[3]): error = 1
                if (win[5] < 0 or win[5] > dim[3]): error = 1
                
                if error == 1:
                    errors += [r, "Range of BC %s is invalid for zone %s."%(n[0],z[0])]
    return errors

def checkBCFaces(t, type):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        dim = Internal.getZoneDim(z)
        r = Internal.getElementRange(z, type="NGON")
        if r is not None: nfaces = r[1]-r[0]+1
        else: nfaces = 0
        nodes = Internal.getNodesFromType2(z, type)
        for n in nodes:
            plist = Internal.getNodesFromName1(n, 'PointList')
            for r in plist:
                faces = r[1]
                faces1 = faces[faces > nfaces]
                faces2 = faces[faces < 1]
                if faces1.size > 0 or faces2.size > 0:
                    errors += [r, "Faces of BC %s is invalid for zone %s."%(n[0],z[0])]
    return errors

#==============================================================================
# Check donor BC ranges
# IN: t: arbre a verifier
# IN: type: type du noeud de BC a verifier (BC_t,...)
# On verifie que le range du donneur est contenu dans la grille donneur
#==============================================================================
def checkDonorRanges(t, type):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType2(z, type)
        for n in nodes:
            donorName = Internal.getValue(n)
            donors = Internal.getNodesFromName2(t, donorName)
            if donors != []:
                dim = Internal.getZoneDim(donors[0])
                r = Internal.getElementRange(donors[0], type="NGON")
                if r is not None: nfaces = r[1]-r[0]+1
                else: nfaces = 0
                prange = Internal.getNodesFromName1(n, 'PointRangeDonor')
                for r in prange:
                    if r[1].shape == (2,3):
                        r[1] = numpy.reshape(r[1], (3,2), order='Fortran')
                    if r[1].shape == (3,2):
                        win = Internal.range2Window(r[1])
                    else: win = [0,0,0,0,0,0] # pas de check en NS
                    error = 0
                    if (win[0] < 0 or win[0] > dim[1]): error = 1
                    if (win[1] < 0 or win[1] > dim[1]): error = 1
                    if (win[2] < 0 or win[2] > dim[2]): error = 1
                    if (win[3] < 0 or win[3] > dim[2]): error = 1
                    if (win[4] < 0 or win[4] > dim[3]): error = 1
                    if (win[5] < 0 or win[5] > dim[3]): error = 1

                    if (error == 1):
                        errors += [r, "Range of donor BC %s is invalid for zone %s."%(n[0],z[0])]
    return errors

def checkDonorFaces(t, type):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType2(z, type)
        for n in nodes:
            donorName = Internal.getValue(n)
            donors = Internal.getNodesFromName2(t, donorName)
            if donors != []:
                plist = Internal.getNodesFromName1(n, 'PointListDonor')
                for r in plist:
                    faces = r[1]
                    faces1 = faces[faces > nfaces]
                    faces2 = faces[faces < 1]
                    if faces1.size > 0 or faces2.size > 0:
                        errors += [r, "Faces of donor BC %s is invalid for zone %s."%(n[0],z[0])]
    return errors

#==============================================================================
# Supprime les BCs de type donne avec des ranges invalides
#==============================================================================
def _correctBCRanges(t, type):
    errors = checkBCRanges(t, type)
    le = len(errors)/2
    for e in xrange(le):
        node = errors[2*e]
        (p, c) = Internal.getParentOfNode(t, node)
        if p is not None:
            (p2, c2) = Internal.getParentOfNode(t, p)
            if p2 is not None: del p2[2][c2]
    return None

#==============================================================================
# Supprime les BCs avec des donor ranges invalides
#==============================================================================
def _correctDonorRanges(t, type):
    errors = checkDonorRanges(t, type)
    le = len(errors)/2
    for e in xrange(le):
        node = errors[2*e]
        (p, c) = Internal.getParentOfNode(t, node)
        if p is not None:
            (p2, c2) = Internal.getParentOfNode(t, p)
            if p2 is not None: del p2[2][c2]
    return None

#==============================================================================
# Verifie les ranges des fenetres opposees
# Le donneur doit exister et les ranges etre coherents
# IN: type: GridConnectivity1to1_t ou GridConnectivity_t
#==============================================================================
def checkOppositRanges(t, type):
    errors = []
    delta = numpy.empty(3, numpy.int32); deltaopp = numpy.empty(3, numpy.int32)
    zones = Internal.getZones(t)
    for z in zones:
        zname = z[0]
        nodes = Internal.getNodesFromType2(z, type)
        dimZone = Internal.getZoneDim(z)[4]
        for n in nodes:
            prange = Internal.getNodesFromName1(n, 'PointRange')
            prangedonor = Internal.getNodesFromName1(n, 'PointRangeDonor')
            mtype = Internal.getNodeFromName1(n, 'GridConnectivityType')
            if mtype is not None:
                mtype = mtype[1]
                if isinstance(mtype, numpy.ndarray): mtype = mtype.tostring()
            else: mtype = 'Match'
            zdonorname = Internal.getValue(n)
            zdonor = Internal.getNodesFromName2(t, zdonorname)
            if zdonor == []:
                errors += [n, "Donor zone %s of BC %s (zone %s) does not exist."%(zdonorname,n[0],z[0])]
            else:
                if mtype != 'Overset' and prange != []:
                    for z in zdonor:
                        if (z[3] == 'Zone_t'): zdonor = z; break
                    # Verifie que le donneur est dans la meme base (trop cher)
                    #(b1, c1) = Internal.getParentOfNode(t, z)
                    #(b2, c2) = Internal.getParentOfNode(t, zdonor)
                    #if (b1[0] != b2[0] and mtype == 'Match'):
                    #    errors += [n, "Donor zone %s of BC %s (zone %s) is not in the same base as zone %s."%(zdonorname,n[0],z[0],z[0])]
                    #else:
                    # Verifie que la paire (prange,prangedonor) existe bien ds la zone zdonor
                    nodesopp = Internal.getNodesFromType2(zdonor, type)
                    dim = Internal.getZoneDim(zdonor)
                    error = 1
                    for nopp in nodesopp:
                        if n is not nopp:  
                            prangeopp = Internal.getNodesFromName1(nopp, 'PointRange')
                            prangedonoropp = Internal.getNodesFromName1(nopp, 'PointRangeDonor')
                            mtypeopp = Internal.getNodesFromName1(nopp, 'GridConnectivityType')
                            if (len(mtypeopp) > 0):
                                mtypeopp = mtypeopp[0]
                            if isinstance(mtypeopp, numpy.ndarray): mtypeopp = mtypeopp.tostring()
                            else: mtypeopp = 'Match'
                            zoppdonorname = nopp[1]
                            if isinstance(zoppdonorname, numpy.ndarray): 
                                zoppdonorname = zoppdonorname.tostring()
                            if zoppdonorname == zname and mtype == mtypeopp:
                                # current zone
                                rangez = Internal.range2Window(prange[0][1]) 
                                rangezd = Internal.range2Window(prangedonor[0][1])
                                # donor zone
                                rangezopp = Internal.range2Window(prangeopp[0][1]) 
                                rangezoppd = Internal.range2Window(prangedonoropp[0][1]) 
                                if rangez == rangezoppd and rangezd == rangezopp: error = 0
                    if (error == 1):
                        errors += [n, "Opposite window from zone %s of BC %s (zone %s) does not exist."%(zdonorname,n[0],zname)]
                    # Check des ranges
                    for ropp in prangedonor:
                        if (ropp[1].shape == (2,3)):
                            ropp[1] = numpy.reshape(ropp[1], (3,2), order='Fortran')
                        if (ropp[1].shape == (3,2)):
                            winopp = Internal.range2Window(ropp[1])
                        else: winopp = [-1,0,0,0,0,0]
                        error = 0
                        if (winopp[0] < 0 or winopp[0] > dim[1]): error = 1
                        if (winopp[1] < 0 or winopp[1] > dim[1]): error = 1
                        if (winopp[2] < 0 or winopp[2] > dim[2]): error = 1
                        if (winopp[3] < 0 or winopp[3] > dim[2]): error = 1
                        if (winopp[4] < 0 or winopp[4] > dim[3]): error = 1
                        if (winopp[5] < 0 or winopp[5] > dim[3]): error = 1
                        if (error == 1):
                            errors += [ropp, "Donor range of BC %s is invalid for zone %s."%(n[0],z[0])]
                        else:
                            if type == 'GridConnectivity1to1_t': # BCMatch only
                                # check consistency of current and donor windows in each direction, taking into account Transform
                                transform = Internal.getNodesFromName1(n, 'Transform')
                                for r in prange:
                                    if (r[1].shape == (2,3)):
                                        r[1] = numpy.reshape(r[1], (3,2), order='Fortran')
                                    if (r[1].shape == (3,2)):
                                        win = Internal.range2Window(r[1])
                                    else: win = [0,0,0,0,0,0]
                                    if transform != []: # not mandatory, [+1,+2,+3] by default.
                                        transform = transform[0][1]
                                    else:
                                        transform = numpy.empty(3, numpy.int32)
                                        transform[0] = 1; transform[1] = 2; transform[2] = 3
                                    delta[0] = abs(win[1] - win[0]) # delta i for win
                                    if dimZone > 1: delta[1] = abs(win[3] - win[2]) # delta j for win
                                    if dimZone == 3: delta[2] = abs(win[5] - win[4]) # delta k for win
                                    deltaopp[0] = abs(winopp[1] - winopp[0]) # delta i for winopp
                                    if dimZone > 1: deltaopp[1] = abs(winopp[3] - winopp[2]) # delta j for winopp
                                    if dimZone == 3: deltaopp[2] = abs(winopp[5] - winopp[4]) # delta k for winopp
                                    if dimZone == 3:
                                        if ((delta[0] != deltaopp[abs(transform[0])-1]) or 
                                            (delta[1] != deltaopp[abs(transform[1])-1]) or
                                            (delta[2] != deltaopp[abs(transform[2])-1])):
                                            errors += [r, "window of BC %s for zone %s does not match with its opposite window."%(n[0],z[0])]
                                    elif dimZone == 2:
                                        if ((delta[0] != deltaopp[abs(transform[0])-1]) or 
                                            (delta[1] != deltaopp[abs(transform[1])-1])):
                                            errors += [r, "window of BC %s for zone %s does not match with its opposite window."%(n[0],z[0])]

                                    elif dimZone == 1:
                                        if delta[0] != deltaopp[abs(transform[0])-1]:
                                                errors += [r, "window of BC %s for zone %s does not match with its opposite window."%(n[0],z[0])]
    return errors

#==============================================================================
# Efface les GCs qui n'ont pas de donneur existant 
#                qui ont des ranges non coherents
#                dont le donneur n'a pas le noeud reciproque
#==============================================================================
def _correctOppositRanges(t, type):
    errors = checkOppositRanges(t, type)
    le = len(errors)/2
    for e in xrange(le):
        node = errors[2*e]
        (p, c) = Internal.getParentOfNode(t, node)
        if p is not None:
            (p2, c2) = Internal.getParentOfNode(t, p)
            if p2 is not None: del p[2][c]
    return None

#==============================================================================
# Si une zone est taggee avec une famille, la famille doit etre declaree
# dans sa base
#==============================================================================
def checkZoneFamily(t):
    zones = Internal.getZones(t)
    errors = []
    for z in zones:
        f = Internal.getNodesFromType1(z, 'FamilyName_t')
        if f != []: # zone taggee
            name = f[0][1]
            if isinstance(name, numpy.ndarray): name = name.tostring()
            # Check for family
            (p, c) = Internal.getParentOfNode(t, z)
            if p is not None:
                ref = Internal.getNodesFromName1(p, name)
                if (ref == []):
                    errors += [p, 'FamilyZone %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
                elif (ref[0][3] != 'Family_t'):
                    errors += [p, 'FamilyZone %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
    return errors

#==============================================================================
# Si une famille n'est pas declaree, on l'ajoute dans la base
#==============================================================================
def _correctZoneFamily(t):
    errors = checkZoneFamily(t)
    le = len(errors)/2
    for e in xrange(le):
        node = errors[2*e]
        name = errors[2*e+1]
        name = name.split(' '); name = name[1]
        (p, c) = Internal.getParentOfNode(t, node)
        if p is not None: C._addFamily2Base(node, name)
    return None

#==============================================================================
# Si une zone utilise une familyBC, la famille doit etre declaree
# dans sa base
#==============================================================================
def checkBCFamily(t):
    zones = Internal.getZones(t)
    errors = []
    for z in zones:
        BCs = Internal.getNodesFromType2(z, 'BC_t')
        for b in BCs:
            f = Internal.getNodesFromType1(b, 'FamilyName_t')
            if f != []: # zone avec BC family
                name = f[0][1]
                if isinstance(name, numpy.ndarray): name = name.tostring()
                # Check for family
                (p, c) = Internal.getParentOfNode(t, z)
                if p is not None:
                    ref = Internal.getNodesFromName1(p, name)
                    if (ref == []):
                        errors += [p, 'FamilyBC %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
                    elif (ref[0][3] != 'Family_t'):
                        errors += [p, 'FamilyBC %s (referenced by zone %s) is not defined in base.'%(name,z[0])]
    return errors

#==============================================================================
# Si une famille BC n'est pas declaree, on l'ajoute dans la base
#==============================================================================
def _correctBCFamily(t):
    errors = checkBCFamily(t)
    le = len(errors)/2
    for e in xrange(le):
        node = errors[2*e]
        name = errors[2*e+1]
        name = name.split(' '); name = name[1]
        (p, c) = Internal.getParentOfNode(t, node)
        if p is not None:
            C._addFamily2Base(node, name, bndType='UserDefined')
    return None

#==============================================================================
# Verifie que dans une base, toutes les zones ont le meme cellDim que la base
#==============================================================================
def checkBaseZonesDim(t):
    errors = []
    bases = Internal.getBases(t)
    for b in bases:
        dimBase = b[1][0]
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            dim = Internal.getZoneDim(z)
            if (dim[4] != dimBase): errors+=[b, "Zone %s cellDim is inconsistent with base %s cellDim."%(z[0],b[0])] 
    return errors

#==============================================================================
# Correct: update la dim de la base si toutes les zones ont le meme
# cellDim, sinon, split la base en plusieurs bases avec des zones
# de meme cellDim
# Ceci est important pour la CGNS lib.
#==============================================================================
def _correctBaseZonesDim(t):
    bases = Internal.getBases(t)
    for b in bases:
        dimBase = b[1][0]
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        z1 = []; z2 = []; z3 = [] # zones de dim 1,2,3                   
        for z in zones:
            dim = Internal.getZoneDim(z)
            if (dim[4] <= 1): z1.append(z)
            elif (dim[4] == 2): z2.append(z)
            elif (dim[4] == 3): z3.append(z)
        lz1 = len(z1); lz2 = len(z2); lz3 = len(z3)
        if (lz1 == 0 and lz2 == 0): Internal.setValue(b, 3)
        elif (lz1 == 0 and lz3 == 0): Internal.setValue(b, 2)
        elif (lz2 == 0 and lz3 == 0): Internal.setValue(b, 1)
        else: print "Bases must be split."
    
    return None

#===============================================================================
# check if the PointRanges for a zone z (ZoneBC ou ZoneGridConnectivity) are 
# compatible with multigrid
#===============================================================================
def checkMGForBCRanges(z,type,multigrid,sizemin):
    puiss = 2**(multigrid)
    errors = []
    nodes = Internal.getNodesFromType2(z, type)
    for n in nodes:
        PRS = Internal.getNodesFromName1(n, 'PointRange')
        for PR in PRS:
            [imin,imax,jmin,jmax,kmin,kmax] = Internal.range2Window(PR[1])
            if imin != imax:
                di = (imax-imin)
                res = (di+1)/puiss
                if di%puiss != 0 : errors+=[n,"BC %s of zone %s is not multigrid of level %d in direction i."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRange for BC %s of zone %s: not enough points on coarse grid in direction i."%(n[0],z[0])]
            if jmin != jmax:
                dj = (jmax-jmin)
                res = (dj+1)/puiss
                if dj%puiss != 0: errors+=[n,"BC %s of zone %s is not multigrid of level %d in direction j."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRange for BC %s of zone %s: not enough points on coarse grid in direction j."%(n[0],z[0])]
            if kmin != kmax:
                dk = (kmax-kmin)
                res = (dk+1)/puiss
                if dk%puiss != 0 : errors+=[n,"BC %s of zone %s is not multigrid of level %d in direction k."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRange for BC %s of zone %s: not enough points on coarse grid in direction k."%(n[0],z[0])]
    return errors

#===============================================================================
# check if the PointRangeDonor for a zone z (ZoneBC ou ZoneGridConnectivity) 
# are compatible with multigrid
#===============================================================================
def checkMGForDonorBCRanges(z, type, multigrid, sizemin):
    puiss = 2**(multigrid)
    errors = []
    nodes = Internal.getNodesFromType2(z, type)
    for n in nodes:
        PRS = Internal.getNodesFromName1(n, 'PointRangeDonor')
        for PR in PRS:
            [imin,imax,jmin,jmax,kmin,kmax] = Internal.range2Window(PR[1])
            if imin != imax:
                di = (imax-imin)
                res = (di+1)/puiss
                if di%puiss != 0: 
                    errors+=[n,"PointRangeDonor for GC %s of zone %s is not multigrid of level %d in direction i."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s: not enough points on coarse grid in direction i."%(n[0],z[0])]
            if jmin != jmax:
                dj = (jmax-jmin)
                res = (dj+1)/puiss
                if dj%puiss != 0: 
                    errors+=[n,"PointRangeDonor for GC %s of zone %s is not multigrid of level %d in direction j."%(n[0],z[0],multigrid)]
                elif res < sizemin:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s: not enough points on coarse grid in direction j."%(n[0],z[0])]
            if kmin != kmax:
                dk = (kmax-kmin)
                res = (dk+1)/puiss
                if dk%puiss != 0:
                    errors+=[n,"PointRangeDonor for GC %s of zone %s is not multigrid of level %d in direction k."%(n[0],z[0],multigrid)]
                elif res < sizemin: errors+=[n,"PointRangeDonor for GC %s of zone %s: not enough points on coarse grid in direction i."%(n[0],z[0])]
    return errors

#==============================================================================
# check if the tree is compatible with multigrid (zones, BC and connectivities)
# level is the MG level that must be ensured: N = 2^level+1
#==============================================================================
def checkMultigrid(t, level=1, nbMinCoarseB=5, nbMinCoarseW=3):
    errors = []
    if level == 0: return errors
    puiss = 2**(level)

    for z in Internal.getZones(t):
        # check Dims
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Structured':
            ni = dims[1]; nj = dims[2]; nk = dims[3]
            res = (ni-1)/puiss
            if (ni-1)%puiss != 0: 
                errors+=[z,"Zone %s is not multigrid of level %d in direction i."%(z[0],level)]
            elif res < nbMinCoarseB:
                errors+=[z,"Zone %s: not enough points on coarse grid for level %d in direction i."%(z[0],level)]
            res = (nj-1)/puiss
            if (nj-1)%puiss != 0: errors+=[z,"Zone %s is not multigrid of level %d in direction j."%(z[0],level)]
            elif res < nbMinCoarseB:
                errors+=[z,"Zone %s: not enough points on coarse grid for level %d in direction j."%(z[0],level)]
            res = (nk-1)/puiss
            if (nk-1)%puiss != 0: errors+=[z,"Zone %s is not multigrid of level %d in direction k."%(z[0],level)]
            elif res < nbMinCoarseB:
                errors+=[z,"Zone %s: not enough points on coarse grid for level %d in direction k."%(z[0],level)]

            # check BC ranges (receptors)
            errors += checkMGForBCRanges(z,'BC_t',level,nbMinCoarseW)
            errors += checkMGForBCRanges(z,'GridConnectivity1to1_t',level,nbMinCoarseW)
            errors += checkMGForBCRanges(z,'GridConnectivity_t',level,nbMinCoarseW)            
            # check BC ranges (donors)
            errors += checkMGForDonorBCRanges(z,'GridConnectivity1to1_t',level,nbMinCoarseW)
            errors += checkMGForDonorBCRanges(z,'GridConnectivity_t',level,nbMinCoarseW)      
    return errors

#=============================================================================
# Check if the number of points of a zone does not exceed sizeMax
#=============================================================================
def checkSize(t, sizeMax=100000000):
    errors = []
    for z in Internal.getZones(t):
        dims = Internal.getZoneDim(z)
        if dims[0] == 'Structured':
            npts = dims[1]*dims[2]*dims[3]
        else: npts = dims[1]
        if npts > sizeMax: errors += [z, "Zone %s exceeds the maximum number of points (Npts=%d)."%(z[0],npts)]
    return errors

#==============================================================================
# Verifie que le type des noeuds est dans CGNSTypes
#==============================================================================
def checkCGNSType(node):
    errors = []
    isStd = Internal.isStdNode(node)
    if (isStd >= 0):
        for c in node[isStd:]: checkCGNSType__(c, errors)
    else: checkCGNSType__(node, errors)
    return errors

def checkCGNSType__(node, errors):
    ntype = node[3]
    if not CGNSTypes.has_key(ntype):
        errors += [node, 'Unknown CGNS type %s for node %s.\n'%(ntype, node[0])]
    sons = node[2]
    for n in sons: checkCGNSType__(n, errors)

#==============================================================================
# Delete les noeuds de type non valide
#==============================================================================
def _correctCGNSType(t):
    errors = checkCGNSType(t)
    le = len(errors)/2
    for e in xrange(le):
        node = errors[2*e]
        (p, c) = Internal.getParentOfNode(t, node)
        if p is not None: del p[2][c]
    return None

#==============================================================================
# Check element nodes  dans t
# Verifie:
# si une zone a NGON+PE et pas de NFace
# si il y a des connectivites multiples
# ou une zone a NGON+PE et pas de NFace
#==============================================================================
def checkElementNodes(t):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        connects = Internal.getElementNodes(z)
        iBE = -1; iBEMultiple = -1; iNGon = -1; iNFace = -1; i = 0
        for c in connects:
            ctype = c[1][0]
            if (ctype == 22): iNGon = i
            elif (ctype == 23): iNFace = i
            else: 
                if (iBE == -1): iBE = i # first
                else: iBEMultiple = 1
            i += 1

        if (iNGon != -1 and iNFace != -1): pass
        elif (iNGon != -1 and iNFace == -1): 
            errors += [z, 'NFace is missing for zone %s.'%z[0]]
        elif (iBEMultiple == 1): 
            errors += [z, 'Multiple BE connectivity for zone %s.'%z[0]]
    return errors

#==============================================================================
# Fait un break connectivity pour les BE multiples
# Ajoute le noeud Face si on a un PE et pas de NFace
#==============================================================================
def _correctElementNodes(t):
    errors = checkElementNodes(t)
    le = len(errors)/2
    for e in xrange(le):
        zone = errors[2*e]
        msg = errors[2*e+1]
        if msg[0:8] == 'Multiple':
            zones = C.breakConnectivity(zone)
            (p,c) = Internal.getParentOfNode(t, zone)
            if p is not None: p[2][c] = zones[0]; p[2] += zones[1:]
        elif msg[0:6] == 'NFace':
            # Look for PE
            PE = Internal.getNodeFromName2(zone, 'ParentElements')
            if PE is not None:
                Internal._adaptPE2NFace(zone, remove=False)
    return None

#===============================================================================
# Check non CGNS varnames in FlowSolution_t and BCDataSet
#===============================================================================
def checkCGNSVarNames(t):
    errors = []
    zones = Internal.getZones(t)
    for z in zones:
        # Change on container
        cont = Internal.getNodesFromType1(z, 'FlowSolution_t')
        for c in cont:
            datas = Internal.getNodesFromType1(c, 'DataArray_t')
            for d in datas:
                n = Internal.getCGNSName(d[0])
                if n != d[0]: errors += [d, '%s not a valid CGNS variable name for zone %s.'%(d[0],z[0])]
        # Change on BCDataSet (if any)
        bcs = Internal.getNodesFromType2(z, 'BC_t')
        for b in bcs:
            datas = Internal.getBCDataSet(z, b)
            for d in datas:
                n = Internal.getCGNSName(d[0])
                if n != d[0]: errors += [d, '%s not a valid CGNS variable name for BCDataSet in zone %s.'%(d[0],z[0])]
    return errors

def _correctCGNSVarNames(t):
    errors = checkCGNSVarNames(t)
    for e in xrange(len(errors)/2):
        node = errors[2*e]
        n = Internal.getCGNSName(node[0])
        node[0] = n
    return None
