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

// Transformations of mesh coordinates ( homothety, contract, symetrize, deform)

# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

extern "C"
{ 
  void k6homothety_(const E_Int& npts,
                    const E_Float* x, const E_Float* y, const E_Float* z,
                    const E_Float& xc, const E_Float& yc, const E_Float& zc,
                    const E_Float& alpha,
                    E_Float* xo, E_Float* yo, E_Float* zo);
  
  void k6deformpoint_(const E_Int& size,
                      const E_Float* dx, const E_Float* dy, const E_Float* dz,
                      const E_Float& amort,
                      const E_Int& ni, const E_Int& nj, const E_Int& nk,
                      E_Float* x, E_Float* y, E_Float* z);
}
// ============================================================================
/* Homothety python from an array describing a mesh */
// ============================================================================
PyObject* K_TRANSFORM::homothety(PyObject* self, PyObject* args)
{
  E_Float xc, yc, zc;
  E_Float alpha;
  PyObject* array;

  if (!PYPARSETUPLEF(args,
                    "O(ddd)d", "O(fff)f",
                    &array, &xc, &yc, &zc, &alpha))
  {
      return NULL;
  }

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;

  res = 
    K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, cn, eltType);

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      if (res == 2) delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "homothety: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
   
    E_Int npts = f->getSize();
    // Homothety
    k6homothety_(npts,
                 f->begin(posx), f->begin(posy), f->begin(posz),
                 xc, yc, zc, alpha,
                 f->begin(posx), f->begin(posy), f->begin(posz));
    
    // Build array
    PyObject* tpl;
    if (res == 1) //structured
    {
      tpl = K_ARRAY::buildArray(*f, varString, nil, njl, nkl);
      delete f;
    } 
    else //unstructured 
    {
      tpl = K_ARRAY::buildArray(*f, varString, *cn, -1, eltType);
      delete f; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "homothety: not a valid array.");
    return NULL;
  }            
}

// ============================================================================
/* Contract python array describing a mesh */
// ============================================================================
PyObject* K_TRANSFORM::contract(PyObject* self, PyObject* args)
{
  E_Float xc, yc, zc;
  E_Float dir1x,dir1y,dir1z;
  E_Float dir2x,dir2y,dir2z;
  E_Float alpha;
  PyObject* array;

  if (!PYPARSETUPLEF(args,
                    "O(ddd)(ddd)(ddd)d", "O(fff)(fff)(fff)f",
                    &array, &xc, &yc, &zc, &dir1x, &dir1y, &dir1z,
                    &dir2x, &dir2y, &dir2z, &alpha))
  {
      return NULL;
  }
  
  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  E_Float norm;

  res = 
    K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, cn, eltType); 

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent( varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent( varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent( varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      if (res == 2) delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "contract: can't find coordinates in array.");
      return NULL;     
    }
    
    posx++; posy++; posz++;

    E_Int npts = f->getSize();

    // Normalisation des vecteurs du plan
    norm = sqrt(dir1x*dir1x + dir1y*dir1y + dir1z*dir1z);
    dir1x = dir1x / norm; 
    dir1y = dir1y / norm; 
    dir1z = dir1z / norm; 
    norm = sqrt(dir2x*dir2x + dir2y*dir2y + dir2z*dir2z);
    dir2x = dir2x / norm; 
    dir2y = dir2y / norm; 
    dir2z = dir2z / norm; 

    // 3 points du plan
    E_Float p1[3];
    E_Float p2[3];
    E_Float p3[3];
    E_Float p[3];
    p1[0] = xc;       p1[1] = yc;       p1[2] = zc;
    p2[0] = xc+dir1x; p2[1] = yc+dir1y; p2[2] = zc+dir1z;
    p3[0] = xc+dir2x; p3[1] = yc+dir2y; p3[2] = zc+dir2z;
    E_Float xint, yint, zint;
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float sigma1, sigma0, dist2;
    E_Boolean in;

    // Contraction
    for (E_Int i = 0; i < npts; i++)
    {
      p[0] = xp[i]; p[1] = yp[i]; p[2] = zp[i];
      K_COMPGEOM::distanceToTriangle(p1, p2, p3, p, 0,
				     dist2, in, xint, yint, zint,
				     sigma0, sigma1);
    
      xp[i] = xint + alpha*(xp[i]-xint); 
      yp[i] = yint + alpha*(yp[i]-yint); 
      zp[i] = zint + alpha*(zp[i]-zint);
    }

    // Build array 
    PyObject* tpl;
    if (res == 1) // structured
    {
      tpl = K_ARRAY::buildArray(*f, varString, 
                                nil, njl, nkl);
      delete f;
    } 
    else // unstructured 
    {
      tpl = K_ARRAY::buildArray(*f, varString,
                                *cn, -1, eltType);
      delete f; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "contract: not a valid array.");
    return NULL;
  }             
}

// ============================================================================
/* Symetry */
/* make a symetric mesh considering the plane passing */
/* by (x0,y0,z0) and of director vectors (v1x,v1y,v1z) */
/* and (v2x,v2y,v2z). */
// ============================================================================
PyObject* K_TRANSFORM::symetrize(PyObject* self, PyObject* args)
{
  E_Float x0, y0, z0;
  E_Float v1x, v1y, v1z;
  E_Float v2x, v2y, v2z;
  PyObject* array;
  
  if (!PYPARSETUPLEF(args,
                    "O(ddd)(ddd)(ddd)", "O(fff)(fff)(fff)",
                    &array, &x0, &y0, &z0, &v1x, &v1y, &v1z, &v2x, &v2y, &v2z))
  {
      return NULL;
  }
  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  E_Float norm;

  res = 
    K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, cn, eltType); 
  
  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent( varString );
    E_Int posy = K_ARRAY::isCoordinateYPresent( varString );
    E_Int posz = K_ARRAY::isCoordinateZPresent( varString );
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      if (res == 2) delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "symetrize: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
   
    E_Int npts = f->getSize();

    // Normalisation des vecteurs du plan
    norm = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
    v1x = v1x / norm; 
    v1y = v1y / norm; 
    v1z = v1z / norm; 
    norm = sqrt(v2x*v2x + v2y*v2y + v2z*v2z);
    v2x = v2x / norm; 
    v2y = v2y / norm;
    v2z = v2z / norm;
   
    // 3 points du plan
    E_Float p1[3];
    E_Float p2[3];
    E_Float p3[3];
    E_Float p[3];
    p1[0] = x0; p1[1] = y0; p1[2] = z0;
    p2[0] = x0+v1x; p2[1] = y0+v1y; p2[2] = z0+v1z;
    p3[0] = x0+v2x; p3[1] = y0+v2y; p3[2] = z0+v2z;
    E_Float xint, yint, zint;
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float sigma1, sigma0, dist2;
    E_Boolean in;

    // Symetry
    for (E_Int i = 0; i < npts; i++)
    {
      p[0] = xp[i]; p[1] = yp[i]; p[2] = zp[i];
      K_COMPGEOM::distanceToTriangle(p1, p2, p3, p, 0,
				     dist2, in, xint, yint, zint,
				     sigma0, sigma1);
    
      xp[i] = 2*xint - xp[i]; 
      yp[i] = 2*yint - yp[i];
      zp[i] = 2*zint - zp[i];
    }

    // Build array 
    PyObject* tpl;
    if (res == 1) //structured
    {
      tpl = K_ARRAY::buildArray(*f, varString, 
                                nil, njl, nkl);
      delete f;
    } 
    else //unstructured 
    {
      tpl = K_ARRAY::buildArray(*f, varString,
                                *cn, -1, eltType);
      delete f; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "symetrize: not a valid array.");
    return NULL;
  }
}
    
// ============================================================================
/* Deform mesh by moving surface of a given vector. */
// ============================================================================
PyObject* K_TRANSFORM::deform(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* normal;

  if (!PyArg_ParseTuple(args, "OO", &array, &normal))
  {
      return NULL;
  }
  
  // Check array
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;

  E_Int res1 = K_ARRAY::getFromArray(array, varString1, f1, 
                                     im1, jm1, km1, cn1, eltType1, true);
  
  E_Int res2 = K_ARRAY::getFromArray(normal, varString2, f2, 
                                     im2, jm2, km2, cn2, eltType2, true);

  // Vecteur et array valides (structure ou non structure) ?
  if (res1 == -1)
  {
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_TypeError,
                    "deform: 1st argument is invalid.");
    return NULL;
  }
  
  if (res2 == -1)
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    PyErr_SetString(PyExc_TypeError,
                    "deform: 2nd argument is invalid.");
    return NULL;
  }

  if (res1 != res2) 
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_TypeError,
                    "deform: 1st and 2nd argument must be both structured or unstructured.");
    return NULL;
  }
  // Presence des coordonnees ?
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_ValueError,
                    "deform: coordinates not found in 1st argument.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  E_Int npts = f1->getSize();
  E_Int nFld = f1->getNfld(); 
  
  // Vecteur normal de dimension 3 ?
  E_Int dim = f2->getNfld();
  if (dim != 3)
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_ValueError,
                    "deform: 2nd argument must have 3 variables defining the vector.");
    return NULL;
  }
  
  // Vecteur et array de la meme taille ?
  if (npts != f2->getSize())
  {
    RELEASESHAREDB(res1, array,f1,cn1);
    RELEASESHAREDB(res2, normal,f2,cn2);
    PyErr_SetString(PyExc_ValueError, 
                    "deform: sizes of 1st and 2nd arguments are not equal.");
    return NULL;
  }

  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  if (res1 == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nFld, varString1, im1, jm1, km1);
  } 
  else //unstructured
  {
    E_Int csize = cn1->getSize()*cn1->getNfld();
    tpl = K_ARRAY::buildArray(nFld, varString1,npts,cn1->getSize(),-1,eltType1,false,csize);      
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn1->begin(), cn1->getSize()*cn1->getNfld());
  }

  E_Float* newFp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF newF(npts,nFld,newFp,true); newF.setAllValuesAt(*f1);// copie

  // Deformation
  E_Float* newFx = newF.begin(posx);
  E_Float* newFy = newF.begin(posy);
  E_Float* newFz = newF.begin(posz);
  E_Float* f1x = f1->begin(posx);
  E_Float* f1y = f1->begin(posy);
  E_Float* f1z = f1->begin(posz);
  E_Float* f2x = f2->begin(1);
  E_Float* f2y = f2->begin(2);
  E_Float* f2z = f2->begin(3);

#pragma omp parallel for default(shared)
  for (E_Int i = 0;  i < npts; i++) 
  {   
    newFx[i] = f1x[i] + f2x[i];
    newFy[i] = f1y[i] + f2y[i];
    newFz[i] = f1z[i] + f2z[i];
  }
  
  RELEASESHAREDB(res1, array,f1,cn1);
  RELEASESHAREDB(res2, normal,f2,cn2);
  return tpl;
} 
// ============================================================================
/* 
   Deform mesh by moving point (xi,yi,zi) of a vector (dx,dy,dz) 
 */
// ============================================================================
PyObject* K_TRANSFORM::deformPoint(PyObject* self, PyObject* args)
{
  E_Float xi, yi, zi;
  E_Float dx, dy, dz;
  E_Float depth;
  E_Float width;
  PyObject* array;

  if (!PYPARSETUPLEF(args,
                    "O(ddd)(ddd)dd", "O(fff)(fff)ff",
                    &array, &xi, &yi, &zi, &dx, &dy, &dz, &depth, &width))
  {
      return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, 
                                    f, im, jm, km, cn, eltType, true);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "deformPoint: invalid array.");
    return NULL;
  }
 
  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "deformPoint: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // normalisation de dx
  E_Float dist;
  dist = dx*dx + dy*dy + dz*dz;
  dist = 1./sqrt(dist);
  dx = dx *dist; dy = dy*dist; dz = dz*dist;
  E_Float sigma = 2./(width*width); 
  
  // Construit l'array resultat
  PyObject* tpl;
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nfld, varString, im, jm, km);
  } 
  else //unstructured 
  {
    E_Int csize = cn->getSize()*cn->getNfld(); 
    tpl = K_ARRAY::buildArray(nfld, varString,
                              npts, cn->getSize(),
                              -1, eltType, false, csize);
  }
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true);
  fn.setAllValuesAt(*f);
  if (res == 2)
  {
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
  }

  // Pointers
  E_Float* xo = fn.begin(posx);
  E_Float* yo = fn.begin(posy);
  E_Float* zo = fn.begin(posz);

  // Regularisation
#pragma omp parallel default(shared)
  {
  E_Float d1, d2, d3, dd, f;
#pragma omp for  
  for (E_Int i = 0; i < npts; i++)
  {
    d1 = xo[i]-xi; d2 = yo[i]-yi; d3 = zo[i]-zi;
    dd = d1*d1+d2*d2+d3*d3;
    f = depth * exp(-dd*sigma);
    xo[i] += dx*f;
    yo[i] += dy*f;
    zo[i] += dz*f;
  }
  }
  
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
