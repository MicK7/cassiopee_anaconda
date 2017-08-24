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
// setBCDataInGhostCellsStruct
# include "converter.h"
using namespace K_FLD;
using namespace std;
//=============================================================================
/* Copy array with ghost values in an array without ghost values */
//=============================================================================
PyObject* K_CONVERTER::setBCDataInGhostCellsStruct(PyObject* self, 
                                                   PyObject* args)
{
  //im,jm,km : dimensions of the zone without ghost cells
  // dataFS : data in zone with ghost cells
  // dataBC : BC data in zone without ghost cells
  // BCRanges : ranges at nodes in real zone (without ghost cells)
  // locI = 0 : solution located at nodes, locI = 1 : solution located at cell centers
  E_Int im, jm, km, d;
  E_Int locI;
  PyObject *dataFS, *dataBC, *BCRanges;// dataFS is dimensioned as ghost cells already im+2*d etc
  if (!PYPARSETUPLEI(args, "OOOlllll", "OOOiiiii", &dataFS, &BCRanges, &dataBC , &im, &jm, &km, &d, &locI)) 
    return NULL;

  // flow field in volume
  FldArrayF* fieldG;
  K_NUMPY::getFromNumpyArray(dataFS, fieldG, true);
  E_Int dim = 3;
  if (km == 1) dim = 2; 
  if (jm == 1) dim = 1;
  //E_Int imjm = im*jm;
  E_Int img = im+2*d;
  E_Int jmg = jm+2*d;
  //E_Int kmg = km+2*d;

  E_Int imgjmg = img*jmg;
  // DataSetInfo : extract list
  E_Int nbc = PyList_Size(BCRanges);
  E_Int nbc0 = PyList_Size(dataBC);
  if ( nbc0 != nbc)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setBCDataInGhostCellsStruct: nb of bc ranges and bc datas are not equal.");
    return NULL;
  }

  vector<FldArrayF*> listOfBCFieldsR;
  FldArrayI rangeBCs(nbc,6);
  for (int v = 0; v < nbc; v++)
  {
    PyObject* tpl = PyList_GetItem(dataBC, v);
    FldArrayF* bcFieldR;
    K_NUMPY::getFromNumpyArray(tpl, bcFieldR, true);    
    listOfBCFieldsR.push_back(bcFieldR);
    PyObject* win = PyList_GetItem(BCRanges, v); // list of integers
    if (PyList_Check(win) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "setBCDataInGhostCellsStruct: each bcrange must be a list.");
      return NULL;
    }
    E_Int sizewin = PyList_Size(win);
    if (sizewin != 6 )
    {
      PyErr_SetString(PyExc_TypeError, 
                      "setBCDataInGhostCellsStruct: window must be [imin,imax,jmin,jmax,kmin,kmax].");
      return NULL;
    }
    E_Int i0 = E_Int(v);
    for (int j = 0; j < sizewin; j++)
    {    
      PyObject* iwin = PyList_GetItem(win, j);
      E_Int j0 = E_Int(j)+1;
      rangeBCs(i0, j0) = PyInt_AsLong(iwin);       
    }
  }

  E_Int* IMIN = rangeBCs.begin(1);
  E_Int* IMAX = rangeBCs.begin(2);
  E_Int* JMIN = rangeBCs.begin(3);
  E_Int* JMAX = rangeBCs.begin(4);
  E_Int* KMIN = rangeBCs.begin(5);
  E_Int* KMAX = rangeBCs.begin(6);

  E_Float* ptrFieldG = fieldG->begin();
  E_Int imin, imax, jmin, jmax, kmin, kmax;  
  if ( locI==0)// nodes : the bcdata directly
  {
    for (E_Int nobc = 0; nobc < nbc; nobc++)
    {
      E_Int indv;
      E_Int noindint=0;
      imin = IMIN[nobc]; jmin = JMIN[nobc]; kmin = KMIN[nobc];
      imax = IMAX[nobc]; jmax = JMAX[nobc]; kmax = KMAX[nobc];
      E_Float* ptrBCFieldR = listOfBCFieldsR[nobc]->begin();

      for (E_Int k = kmin-1; k < kmax; k++)
        for (E_Int j = jmin-1; j < jmax; j++)
          for (E_Int i = imin-1; i < imax; i++)
          {
            indv = (i+d)+(j+d)*img+(k+d)*imgjmg; 
            ptrFieldG[indv] = ptrBCFieldR[noindint];
            noindint++;
          }
    }
  }
  else //centers
  {
    E_Int imgc = img-1;// dim ni at centers for ghost cells zone
    E_Int jmgc = jmg-1;
    //E_Int kmgc = kmg-1;
    E_Int imgcjmgc=imgc*jmgc;
    E_Int incr;
    if ( dim == 3)
    {
      for (E_Int nobc = 0; nobc < nbc; nobc++)
      {
        E_Int indadj, indadj0, indghost, ig, jg, kg;
        E_Int noindint=0;
        imin = IMIN[nobc]+d; jmin = JMIN[nobc]+d; kmin = KMIN[nobc]+d;
        imax = IMAX[nobc]+d; jmax = JMAX[nobc]+d; kmax = KMAX[nobc]+d;
        E_Float* ptrBCFieldR = listOfBCFieldsR[nobc]->begin();
        if ( imin == imax)
        {
          if ( imin == d+1 )
          { 
            incr = -1;
            ig = (imin-1);
          }
          else 
          {
            incr = 1;
            ig = (imin-2);
          }          
          for (E_Int k = kmin-1; k < kmax-1; k++)
            for (E_Int j = jmin-1; j < jmax-1; j++)
            {
              indadj0 = ig+j*imgc+k*imgcjmgc;
              for (E_Int d0 = 0; d0 < d; d0++)
              {
                indadj   = indadj0-incr*d0;
                indghost = indadj0+incr*(d0+1);
                ptrFieldG[indghost] = 2.*ptrBCFieldR[noindint]-ptrFieldG[indadj];
              }
              noindint++;
            }
        }
        else if ( jmin == jmax)
        {
          if (jmin == d+1) 
          { 
            incr = -imgc;            
            jg = (jmin-1);
          }
          else 
          { 
            incr = imgc; 
            jg = (jmin-2); 
          } 
          for (E_Int k = kmin-1; k < kmax-1; k++)
            for (E_Int i = imin-1; i < imax-1; i++)
            {
              indadj0 = i + jg*imgc+k*imgcjmgc;
              for (E_Int d0 = 0; d0 < d; d0++)
              {
                indadj   = indadj0-incr*d0;
                indghost = indadj0+incr*(d0+1);
                ptrFieldG[indghost] = 2.*ptrBCFieldR[noindint]-ptrFieldG[indadj];
              }
              noindint++;
            }
        }
        else if ( kmin == kmax)
        {
          if (kmin == d+1) 
          { 
            incr = -imgcjmgc;            
            kg = (kmin-1);
          }
          else 
          { 
            incr = imgcjmgc; 
            kg = (kmin-2); 
          } 

          for (E_Int j = jmin-1; j < jmax-1; j++)
            for (E_Int i = imin-1; i < imax-1; i++)
            {
              indadj0 = i + j*imgc+kg*imgcjmgc;
              for (E_Int d0 = 0; d0 < d; d0++)
              {
                indadj   = indadj0-incr*d0;
                indghost = indadj0+incr*(d0+1);
                ptrFieldG[indghost] = 2.*ptrBCFieldR[noindint]-ptrFieldG[indadj];
              }
              noindint++;
            }
        }                  
      }
    }//dim3
    else 
    {
      for (E_Int nobc = 0; nobc < nbc; nobc++)
      {
        E_Int indadj, indadj0, indghost, ig, jg;
        E_Int noindint=0;
        imin = IMIN[nobc]+d; jmin = JMIN[nobc]+d;
        imax = IMAX[nobc]+d; jmax = JMAX[nobc]+d; 
        E_Float* ptrBCFieldR = listOfBCFieldsR[nobc]->begin();
        if ( imin == imax)
        {
          if ( imin == d+1 )
          { 
            incr = -1;
            ig = (imin-1);
          }
          else 
          {
            incr = 1;
            ig = (imin-2);
          }
          for (E_Int j = jmin-1; j < jmax-1; j++)
          {
            indadj0 = ig+j*imgc;
            for (E_Int d0 = 0; d0 < d; d0++)
            {
              indadj   = indadj0-incr*d0;
              indghost = indadj0+incr*(d0+1);
              ptrFieldG[indghost] = 2.*ptrBCFieldR[noindint]-ptrFieldG[indadj];
            }
            noindint++;
          }
        }
        else if ( jmin == jmax)
        {
          if (jmin == d+1) 
          { 
            incr = -imgc;            
            jg = (jmin-1);
          }
          else 
          { 
            incr = imgc; 
            jg = (jmin-2); 
          } 

          for (E_Int i = imin-1; i < imax-1; i++)
          {
            indadj0 = i + jg*imgc;

            for (E_Int d0 = 0; d0 < d; d0++)
            {
              indadj   = indadj0-incr*d0;
              indghost = indadj0+incr*(d0+1);
              ptrFieldG[indghost] = 2.*ptrBCFieldR[noindint]-ptrFieldG[indadj];
            }
            noindint++;
          }
        }
      }
    }
  }

  Py_INCREF(Py_None);
  return Py_None;  
}
