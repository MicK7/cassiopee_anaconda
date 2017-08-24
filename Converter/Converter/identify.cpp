/*    
    Copyright 2013-2017 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

// Routines d'identification geometrique

# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Identifie les noeuds de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: zone a identfier
   Retourne la liste des noeuds dans la numerotation du KDT. */
// ============================================================================
PyObject* K_CONVERTER::identifyNodes(PyObject* self, PyObject* args)
{  
  PyObject* array; PyObject* hook; E_Float tol;
  if (!PYPARSETUPLEF(args, "OOd", "OOf", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3 &&
      *type != 100 && *type != 102 && *type != 103)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: array is invalid.");
    return NULL;
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res,array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int npts = f->getSize();
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);

  PyObject* ac = K_NUMPY::buildNumpyArray(npts, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Remplissage
#pragma omp parallel default(shared)
  {
  E_Float pt[3];
  E_Float xf,yf,zf,dx,dy,dz;
  E_Int ind;
#pragma omp for schedule(dynamic)
  for (E_Int i = 0; i < npts; i++)
  {
    xf = xp[i]; yf = yp[i]; zf = zp[i];
    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    ind = globalKdt->getClosest(pt); // closest pt
    dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
    //d = dx*dx + dy*dy + dz*dz;
    //if (d < tol) nptr[i] = ind+1;
    //else nptr[i] = -1;
    if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
    else nptr[i] = -1;
  }
  }
  RELEASESHAREDB(res, array, f, cnl);
  return ac;
}

// ============================================================================
/* Identifie les centres des faces de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne la liste des faces. */
// ============================================================================
PyObject* K_CONVERTER::identifyFaces(PyObject* self, PyObject* args)
{ 
  PyObject* array; PyObject* hook;
  E_Float tol;
  if (!PYPARSETUPLEF(args, "OOd", "OOf", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3&&
      *type != 100 && *type != 102 && *type != 103)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must be a NGON.");
    return NULL; 
  }
  if (strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must be a NGON.");
    return NULL; 
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int nfaces = (*cnl)[0];
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);
  E_Int* cnp = cnl->begin();

  PyObject* ac = K_NUMPY::buildNumpyArray(nfaces, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Remplissage
  FldArrayI posFace;
  K_CONNECT::getPosFaces(*cnl, posFace);
  E_Int* posFacep = posFace.begin();

#pragma omp parallel default(shared)
  {
  E_Int nv, ind, pos;
  E_Float xf, yf, zf, inv, dx, dy, dz;
  E_Float pt[3];
  E_Int* ptr;
  
#pragma omp for schedule(dynamic)
  for (E_Int i = 0; i < nfaces; i++)
  {
    pos = posFacep[i];
    ptr = &cnp[pos];
    nv = ptr[0];
    xf = 0.; yf = 0.; zf = 0.;

    for (E_Int n = 1; n <= nv; n++)
    {
      ind = ptr[n]-1;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
    }
    inv = 1./nv; xf *= inv; yf *= inv; zf *= inv;

    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    ind = globalKdt->getClosest(pt); // closest pt
    dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
    //d = dx*dx + dy*dy + dz*dz;
    //if (d < tol) nptr[i] = ind+1; 
    //else nptr[i] = -1;
    if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1; 
    else nptr[i] = -1;
  }
  }
  RELEASESHAREDU(array, f, cnl);
  return ac;
}

// ============================================================================
/* Identifie les centres des elements de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne les indices des points du kdtree correspondant. */
// ============================================================================
PyObject* K_CONVERTER::identifyElements(PyObject* self, PyObject* args)
{ 
  PyObject* array; PyObject* hook;
  E_Float tol;
  if (!PyArg_ParseTuple(args, "OOd", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3&&
      *type != 100 && *type != 102 && *type != 103) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: array must be a NGON.");
    return NULL; 
  }
  
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;
  PyObject* ac = NULL;

  if (strcmp(eltType, "NGON") == 0)
  {  
    // Cree le numpy de sortie
    E_Int sizef = (*cnl)[1];
    E_Int nelts = (*cnl)[sizef+2];
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);

    ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

    // Remplissage
    E_Int* ptr = cnl->begin();
    FldArrayI posFace;
    K_CONNECT::getPosFaces(*cnl, posFace);
    E_Int* posFacep = posFace.begin();
    FldArrayI posElts;
    K_CONNECT::getPosElts(*cnl, posElts);
    E_Int* posEltsp = posElts.begin();

    //printf("nelts=%d, tol=%g\n", nelts, tol);

#pragma omp parallel default(shared)
    {
    E_Float xf,yf,zf,dx,dy,dz;
    E_Int ind, indp, nf, pos, nv, c;
    E_Float inv;
    E_Float pt[3];
    E_Int* ptrFace;
    E_Int* ptrElt;
    
#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < nelts; i++)
    {
      ptrElt = &ptr[posEltsp[i]];
      nf = ptrElt[0];
      xf = 0.; yf = 0.; zf = 0.; c = 0;

      for (E_Int n = 1; n <= nf; n++)
      { 
        ind = ptrElt[n]-1;
        pos = posFacep[ind];
        ptrFace = &ptr[pos];
        nv = ptrFace[0];
        for (E_Int p = 1; p <= nv; p++)
        {
          indp = ptrFace[p]-1; xf += xp[indp]; yf += yp[indp]; zf += zp[indp]; c++;
        }
      }
      inv = 1./E_Float(c); xf *= inv; yf *= inv; zf *= inv;

      pt[0] = xf; pt[1] = yf; pt[2] = zf;
      ind = globalKdt->getClosest(pt); // closest pt
      dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
      //d = (xt[ind]-xf)*(xt[ind]-xf)+(yt[ind]-yf)*(yt[ind]-yf)+(zt[ind]-zf)*(zt[ind]-zf);
      //if (d < tol) nptr[i] = ind+1;
      //else nptr[i] = -1;
      if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
      else nptr[i] = -1;
      //if (K_FUNC::E_abs(xf-3.)<1.e-1 && K_FUNC::E_abs(yf-9.5)<1.e-1 && K_FUNC::E_abs(zf-4.5)<1.e-1)
      //if (nptr[i] != -1)
      //printf("HERE %d: %f %f %f -> %f %f %f -> %d\n",i,xf,yf,zf,dx,dy,dz,nptr[i]);
    }
    }
  }
  else if (strcmp(eltType,"BAR")==0 || strcmp(eltType,"TRI")==0 || 
           strcmp(eltType,"QUAD")==0 || strcmp(eltType,"TETRA")==0 || 
           strcmp(eltType,"HEXA")==0 || strcmp(eltType,"PENTA")==0 )
  {
    E_Int nelts = cnl->getSize();
    E_Int nvert = cnl->getNfld();
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);
    ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);
    E_Float pt[3];
    E_Float inv = 1./nvert;

#pragma omp parallel for default(shared) private(pt) schedule(dynamic)
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Float xf,yf,zf,dx,dy,dz;
      E_Int ind;
      xf = 0.; yf = 0.; zf = 0.; 
      for (E_Int n = 1; n <= nvert; n++)
      {
        ind = (*cnl)(i,n);
        xf += xp[ind]; yf += yp[ind]; zf += zp[ind];        
      }
      pt[0] = xf*inv; pt[1] = yf*inv; pt[2] = zf*inv;
      ind = globalKdt->getClosest(pt); // closest pt
      dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
      //d = dx*dx + dy*dy + dz*dz;
      //if (d < tol) nptr[i] = ind+1;
      //else nptr[i] = -1;
      if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
      else nptr[i] = -1;
    }
  }
  else
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: invalid type of array.");
    return NULL; 
  }
  RELEASESHAREDU(array, f, cnl);
  return ac;
}
