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
// compute variables such as p, primitive variables

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include "post.h"

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6computevelocity_(const E_Int& npts, const E_Int& neq,
			  const E_Int& posro,
			  const E_Int& posrou, const E_Int& posrov,
			  const E_Int& posrow,
			  const E_Float* field, 
			  E_Float* velo);

  void k6computepressure_(const E_Int& cellt,const E_Float& gamma, 
                          const E_Float* velo, 
                          const E_Float* ro, const E_Float* roe,
                          E_Float* p);
}

//=============================================================================
/* Compute variables */
//=============================================================================
PyObject* K_POST::computeVariables(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars0;
  E_Float gamma, rgp, s0, betas, Cs, mus, Ts;
  if (!PYPARSETUPLEF(args,
                    "OOddddddd", "OOfffffff",
                    &array, &vars0, &gamma, &rgp, &s0, &betas, &Cs, &mus, &Ts))
  {
      return NULL;
  }

  // calcul de betas donne par mus,Ts
  if (K_FUNC::fEqualZero(mus) == false && 
      K_FUNC::fEqualZero(Ts) == false)
  {
    betas = mus*(Ts+Cs)/(Ts*sqrt(Ts));
    //printf("Info: computeVariables: betas = %12.15e\n", betas);
  }

  // Check array
  char* varString; char* eltType;
  char varStringOut[K_ARRAY::VARSTRINGLENGTH]; varStringOut[0] = '\0';
  FldArrayF* f; FldArrayI* c;
  E_Int res; 
  E_Int ni, nj, nk; // number of points of array
  res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, c, 
                              eltType, true);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeVariables: unknown type of array.");
    return NULL;
  } 
     
  // Extrait les variables a calculer de la chaine vars0. 
  // Insert dans vars uniquement celles qui seront effectivement calculees
  vector<char*> vars;
  if (PyList_Check(vars0) != 0)
  {
    for (int i = 0; i < PyList_Size(vars0); i++)
    {
      PyObject* tpl0 = PyList_GetItem(vars0, i);
      if (PyString_Check(tpl0) == 0)
      {
        printf("Warning: computeVariables: varname must be a string. Skipped...\n");
      }
      else 
      {
        char* str = PyString_AsString(tpl0);
        char tmpVarString[K_ARRAY::VARSTRINGLENGTH];
        short ok = checkAndExtractVariables(str, vars, tmpVarString);
        if (ok != 0) 
        {
          if (varStringOut[0] == '\0') strcpy(varStringOut, tmpVarString);
          else {strcat(varStringOut, ","); strcat(varStringOut, tmpVarString);}
        }
      }
    }
  }
  else 
  {
    if (PyString_Check(vars0) == 0) printf("Warning: computeVariables: varname must be a string. Skipped...\n");
    else 
    {  
      char* str = PyString_AsString(vars0);
      checkAndExtractVariables(str, vars, varStringOut);
    }
  }
  E_Int nvarout = vars.size();// variables a calculer
  if (nvarout == 0)
  {
    RELEASESHAREDB(res, array, f, c); 
    PyErr_SetString(PyExc_TypeError,
                    "computeVariables: no variable computed.");
    return NULL;
  }

  // Les variables conservatives sont requises
  E_Int posro, posrou, posrov, posrow, posroe;
  posro = K_ARRAY::isDensityPresent(varString); 
  posrou = K_ARRAY::isMomentumXPresent(varString); 
  posrov = K_ARRAY::isMomentumYPresent(varString); 
  posrow = K_ARRAY::isMomentumZPresent(varString); 
  posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);   
  if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
  {
    printf("Warning: computeVariables: one conservative field was not found. Variables are not computed.\n");
    RELEASESHAREDB(res, array, f, c); return array;
  }
  posro++; posrou++; posrov++; posrow++; posroe++;

  FldArrayF* fnew = new FldArrayF(f->getSize(), nvarout);

  // calcul des variables composees
  E_Int ok = computeCompVariables(
    *f, posro, posrou, posrov, posrow, posroe, 
    gamma, rgp, s0, betas, Cs, vars, *fnew); 

  // nettoyage...
  E_Int varsSize = vars.size();
  for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];

  if (ok == 0) //erreur de developpt
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeVariables: invalid string.\n");
    delete fnew; RELEASESHAREDB(res, array, f, c); return NULL;
  }

  if (res == 1)
  {
    PyObject* tpl = K_ARRAY::buildArray(*fnew, varStringOut, 
                                        ni, nj, nk);
    delete fnew; 
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else
  {
    PyObject* tpl = K_ARRAY::buildArray(*fnew, varStringOut, 
                                        *c, -1, eltType);
    delete fnew; 
    RELEASESHAREDU(array, f, c);
    return tpl;
  }
}
//-----------------------------------------------------------------------------
/* Extrait de la chaine vars0 les variables a calculer. Une verification 
   est effectuee sur les noms de variables. La chaine varStringOut est 
   aussi construite pour l'array de sortie 
   IN: vars0: chaine contenant les variables a extraire
   OUT: vars: vecteur contenant les variables a calculer
   OUT: varStringOut: chaine de variables calculees pour l'array de sortie 
   retourne 0 si aucune variable n a ete trouvee.
*/
//-----------------------------------------------------------------------------
short K_POST::checkAndExtractVariables(char* vars0, vector<char*>& vars, 
                                       char* varStringOut)
{
  vector<char*> tmpvars;
  //extraction de toutes les variables
  K_ARRAY::extractVars(vars0, tmpvars);

  E_Int maxStringLength = 80;

  E_Int first = 1; // test si premier element de vars
  // verification
  E_Int tmpvarsSize = tmpvars.size();
  short ok = 0;

  for (E_Int v = 0 ; v < tmpvarsSize; v++)
  {
    if (K_STRING::cmp(tmpvars[v], "VelocityX") == 0 ||
        K_STRING::cmp(tmpvars[v], "VelocityY") == 0 ||
        K_STRING::cmp(tmpvars[v], "VelocityZ") == 0 ||
        K_STRING::cmp(tmpvars[v], "VelocityMagnitude") == 0 ||
        K_STRING::cmp(tmpvars[v], "Pressure") == 0 ||
        K_STRING::cmp(tmpvars[v], "Temperature") == 0 ||
        K_STRING::cmp(tmpvars[v], "Entropy") == 0 ||
        K_STRING::cmp(tmpvars[v], "Enthalpy") == 0 ||
        K_STRING::cmp(tmpvars[v], "Mach") == 0 ||
        K_STRING::cmp(tmpvars[v], "ViscosityMolecular") == 0 ||
        K_STRING::cmp(tmpvars[v], "PressureStagnation") == 0  ||
        K_STRING::cmp(tmpvars[v], "TemperatureStagnation") == 0 ||
        K_STRING::cmp(tmpvars[v], "PressureDynamic") == 0 
//         K_STRING::cmp(tmpvars[v],"RotatingVelocityX") == 0 ||
//         K_STRING::cmp(tmpvars[v],"RotatingVelocityY") == 0 ||
//         K_STRING::cmp(tmpvars[v],"RotatingVelocityZ") == 0 ||
//         strcmp(tmpvars[v],"RotatingVelocityMagnitude") == 0
      )
    {
      char* tmp = new char[maxStringLength];
      strcpy(tmp, tmpvars[v]);
      vars.push_back(tmp);
      if (first == 1)
      {
        strcpy(varStringOut, tmpvars[v]);
        first = 0;
        ok = 1;
      }
      else 
      {
        strcat(varStringOut, ",");
        strcat(varStringOut, tmpvars[v]);
      }
    }
    else 
    {
      printf("Warning: computeVariables: %s is not a valid string name. Skipped...\n", tmpvars[v]);
    }
  }
  for (E_Int v = 0; v < tmpvarsSize; v++) delete [] tmpvars[v];
  
  return ok;
}

//-----------------------------------------------------------------------------
// Calcule les variables composees
//-----------------------------------------------------------------------------
E_Int K_POST::computeCompVariables(const FldArrayF& f, const E_Int posro,
                                   const E_Int posrou, const E_Int posrov,
                                   const E_Int posrow, const E_Int posroe,
                                   const E_Float gamma, const E_Float rgp,
                                   const E_Float s0,  
                                   const E_Float betas, const E_Float Cs,
                                   vector<char*>& vars,
                                   FldArrayF& fnew)
{ 
  E_Int npts = f.getSize();

  FldArrayF velo;// vitesse absolue 
  FldArrayF press;// pression statique
  FldArrayF temp; //temperature statique
  FldArrayF mu; //viscosite du fluide
  FldArrayF mach;// mach absolu

  E_Int varsSize = vars.size();
  for (E_Int v = 0 ; v < varsSize; v++)
  {
    E_Float* fnewv = fnew.begin(v+1);
    if (K_STRING::cmp(vars[v], "VelocityX") == 0) //vitesse absolue vx 
    {
      computeVelocity(f, posro, posrou, posrov, posrow, velo);
      E_Float* velox = velo.begin(1);
      for (E_Int ind = 0; ind < npts; ind++) fnewv[ind] = velox[ind];
    }
    else if (K_STRING::cmp(vars[v], "VelocityY") == 0) //vitesse absolue vy
    {
      computeVelocity(f, posro, posrou, posrov, posrow, velo);
      E_Float* veloy = velo.begin(2);
      for (E_Int ind = 0; ind < npts; ind++) fnewv[ind] = veloy[ind];
    }
    else if (K_STRING::cmp(vars[v], "VelocityZ") == 0) //vitesse absolue vz
    {
      computeVelocity(f, posro, posrou, posrov, posrow, velo);
      E_Float* veloz = velo.begin(3);
      for (E_Int ind = 0; ind < npts; ind++) fnewv[ind] = veloz[ind];
    }
    else if (K_STRING::cmp(vars[v], "VelocityMagnitude") == 0)// vitesse abs: module
    {
      computeVelocity(f, posro, posrou, posrov, posrow, velo);
      E_Float* velox = velo.begin(1);
      E_Float* veloy = velo.begin(2);
      E_Float* veloz = velo.begin(3);
      E_Float u1, u2, u3;
      for (E_Int ind = 0; ind < npts; ind++)
      {
        u1 = velox[ind]; u2 = veloy[ind]; u3 = veloz[ind];
        fnewv[ind] = sqrt(u1*u1+u2*u2+u3*u3);
      }
    }
    else if (K_STRING::cmp(vars[v], "Pressure") == 0) // pression statique 
    {
      computePressure(f, posro, posrou, posrov, posrow, posroe,
                      gamma, velo, press);
      E_Float* pressp = press.begin();
      for (E_Int ind = 0; ind < npts; ind++) fnewv[ind] = pressp[ind];
    }
    else if (K_STRING::cmp(vars[v], "Temperature") == 0) //temperature statique
    {
      computeTemperature(f, posro, posrou, posrov, posrow, posroe,
                         gamma, rgp, velo, press, temp);
      E_Float* tempp = temp.begin();
      for (E_Int ind = 0 ; ind < npts; ind++) fnewv[ind] = tempp[ind];
    }
    else if (K_STRING::cmp(vars[v], "Entropy") == 0)
    {
      FldArrayF s(npts);
      computeEntropy(f, posro, posrou, posrov, posrow, posroe,
                     gamma, rgp, s0, velo, temp, press, s);
      E_Float* sp = s.begin();
      for (E_Int ind = 0 ; ind < npts; ind++) fnewv[ind] = sp[ind];
    }
    else if (K_STRING::cmp(vars[v], "Enthalpy") == 0)
    {
      FldArrayF h(npts);
      computeEnthalpy(f, posro, posrou, posrov, posrow, posroe,
                      gamma, velo, press, h);
      E_Float* hp = h.begin();
      for (E_Int ind = 0 ; ind < npts; ind++) fnewv[ind] = hp[ind];  
    }
    else if (K_STRING::cmp(vars[v], "Mach") == 0)
    {
      computeMach(f, posro, posrou, posrov, posrow, posroe, gamma,
                  velo, press, mach);
      E_Float* machp = mach.begin();
      for (E_Int ind = 0 ; ind < npts; ind++) fnewv[ind] = machp[ind];
    }
    else if (K_STRING::cmp(vars[v], "ViscosityMolecular") == 0) //viscosite du fluide 
    {
      computeMu(f, posro, posrou, posrov, posrow, posroe,
		gamma, rgp, betas, Cs, velo, press, temp, mu);
      E_Float* mup = mu.begin();
      for (E_Int ind = 0; ind < npts; ind++) fnewv[ind] = mup[ind];
    }
    else if (K_STRING::cmp(vars[v], "PressureStagnation") == 0) //pression d'arret
    {
      computePressure(f, posro, posrou, posrov, posrow, posroe,
                      gamma, velo, press);
      computeMach(f, posro, posrou, posrov, posrow, posroe, gamma,
                  velo, press, mach);
      E_Float coef = gamma/(gamma-1.);
      E_Float coef1 = 0.5*(gamma-1);
      E_Float* pressp = press.begin();
      E_Float* machp = mach.begin();
      for (E_Int ind = 0; ind < npts; ind++)
      {
        fnewv[ind] = pressp[ind]*pow(1.+coef1*machp[ind]*machp[ind], coef);
      }
    }
    else if (K_STRING::cmp(vars[v], "TemperatureStagnation") == 0) //temperature d'arret
    {
      computeTemperature(f, posro, posrou, posrov, posrow, posroe,
                         gamma, rgp, velo, press, temp);
      computeMach(f, posro, posrou, posrov, posrow, posroe, gamma,
                  velo, press, mach);
      E_Float coef1 = 0.5*(gamma-1);
      E_Float* tempp = temp.begin();
      E_Float* machp = mach.begin();
      for (E_Int ind = 0; ind < npts; ind++)
        fnewv[ind] = tempp[ind]*(1.+coef1*machp[ind]*machp[ind]);        
    }
    else if (K_STRING::cmp(vars[v], "PressureDynamic") == 0) //pression dynamique 
    {
      computePressure(f, posro, posrou, posrov, posrow, posroe,
                      gamma, velo, press);
      computeMach(f, posro, posrou, posrov, posrow, posroe, gamma,
                  velo, press, mach);
      const E_Float* pressp = press.begin();
      const E_Float* machp = mach.begin();
      // pression dynamique 
      for (E_Int ind = 0; ind < npts; ind++)
        fnewv[ind] = 0.5*gamma*pressp[ind]*machp[ind]*machp[ind];

    }
    else return 0;
  }
  return 1;
}

//=============================================================================
// Calcul de la vitesse absolue
// IN: f: champ (doit contenir les variables conservatives)
// OUT: velo: champ de vitesse
//=============================================================================
void K_POST::computeVelocity(const FldArrayF& f, 
                             const E_Int posro, const E_Int posrou, 
                             const E_Int posrov, const E_Int posrow,
                             FldArrayF& velo)
{
  if (velo.getSize() != 0) return; // deja calcule
 
  E_Int npts = f.getSize();
  E_Int nfld = f.getNfld();

  velo.malloc(npts, 3);
  k6computevelocity_(npts, nfld,
  		     posro, posrou, posrov, posrow,
  		     f.begin(), velo.begin());
}

//=============================================================================
// Calcul de la pression statique
// IN: f: champ (doit contenir les variables conservatives)
// IN: gamma: constante gaz parfait
// OUT: velo: champ de vitesse
// OUT: press: pression statique pour un gaz parfait
//=============================================================================
void K_POST::computePressure(const FldArrayF& f, 
                             const E_Int posro, const E_Int posrou, 
                             const E_Int posrov, const E_Int posrow,
                             const E_Int posroe, const E_Float gamma,
                             FldArrayF& velo, FldArrayF& press)
{
  if (press.getSize() != 0) return; // deja calcule
  
  E_Int npts = f.getSize();
  press.malloc(npts);
    
  computeVelocity(f, posro, posrou, posrov, posrow, velo);
  
  // calcul de p = gam1 * roe, e energie interne
  k6computepressure_(npts, gamma, velo.begin(),
		     f.begin(posro), f.begin(posroe), 
		     press.begin());
}

//=============================================================================
// Calcul de la temperature statique
// IN: f: champ (doit contenir les variables conservatives)
// IN: gamma, rgp
// OUT: velo, press, temp
//=============================================================================
void K_POST::computeTemperature(const FldArrayF& f, 
                                const E_Int posro, const E_Int posrou, 
                                const E_Int posrov, const E_Int posrow,
                                const E_Int posroe, const E_Float gamma,
                                const E_Float rgp,
                                FldArrayF& velo, FldArrayF& press,
                                FldArrayF& temp)
{
  if (temp.getSize() != 0) return; // deja calcule
  
  E_Int npts = f.getSize();
  temp.malloc(npts);
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press); 
  const E_Float* rop = f.begin(posro);
  E_Float* tempp = temp.begin();
  E_Float* pressp = press.begin();
  for (E_Int ind = 0; ind < npts; ind++)
    tempp[ind] = pressp[ind] / (rop[ind]*rgp);
}

//=============================================================================
// Calcul de la viscosite du fluide
//=============================================================================
void K_POST::computeMu(const FldArrayF& f, 
                       const E_Int posro, const E_Int posrou, 
                       const E_Int posrov, const E_Int posrow,
                       const E_Int posroe, 
                       const E_Float gamma, const E_Float rgp,
                       const E_Float betas, const E_Float Cs,
                       FldArrayF& velo, FldArrayF& press,
                       FldArrayF& temp, FldArrayF& mu)
{
  if (mu.getSize() != 0) return; // deja calcule
  
  E_Float t, inv;
  E_Int npts = f.getSize();
  mu.malloc(npts);
  computeTemperature(f, posro, posrou, posrov, posrow, posroe,
                     gamma, rgp, velo, press, temp);
  E_Float* tempp = temp.begin();
  E_Float* mup = mu.begin();
  for (E_Int ind = 0 ; ind < npts; ind++)
  {
    t = tempp[ind];
    inv = 1.+Cs/t;      
    mup[ind] = betas * sqrt(t) / inv;
  }
}

//=============================================================================
/* Calcul de l'entropie: s = s0 + rgp*gam/(gam-1)*log(T) - rgp*log(P)
   ou s0 = sref -  rgp*gam/(gam-1)*log(Tref) + rgp*log(Pref)
   IN: f: champ (doit contenir les variables conservatives)
   IN: gamma, rgp: constantes des gaz parfaits
*/
//=============================================================================
void K_POST::computeEntropy(const FldArrayF& f, 
                            const E_Int posro, const E_Int posrou, 
                            const E_Int posrov, const E_Int posrow,
                            const E_Int posroe, const E_Float gamma,
                            const E_Float rgp, const E_Float s0,
                            FldArrayF& velo,
                            FldArrayF& temp, FldArrayF& press,
                            FldArrayF& s)
{ 
  E_Int npts = f.getSize();
  E_Float cp = rgp * gamma / (gamma-1.);
  
  //calcul de la pression si pas deja calculee
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press); 
  //calcul de la temperature si pas deja calculee
  computeTemperature(f, posro, posrou, posrov, posrow, posroe,
                    gamma, rgp, velo, press, temp);
  E_Float* sp = s.begin();
  E_Float* tempp = temp.begin();
  E_Float* pressp = press.begin();
  for (E_Int ind = 0 ; ind < npts; ind++)
    sp[ind] = s0 + cp * log(tempp[ind]) - rgp * log(pressp[ind]);
}
//=============================================================================
// calcul de l'enthalpie
//=============================================================================
void K_POST::computeEnthalpy(const FldArrayF& f, 
                             const E_Int posro, const E_Int posrou, 
                             const E_Int posrov, const E_Int posrow,
                             const E_Int posroe, const E_Float gamma,
                             FldArrayF& velo, FldArrayF& press, FldArrayF& h)
{
  E_Int npts = f.getSize();

  // calcul de la pression si pas deja calculee
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press);
  E_Float gam6 = gamma/(gamma-1.);
  const E_Float* rop = f.begin(posro);
  E_Float* hp = h.begin();
  E_Float* pressp = press.begin();

  for (E_Int i = 0; i < npts; i++) hp[i] = gam6 * pressp[i]/rop[i];
}

//=============================================================================
// Calcul du mach 
//=============================================================================
void K_POST::computeMach(const FldArrayF& f, 
                         const E_Int posro, const E_Int posrou, 
                         const E_Int posrov, const E_Int posrow,
                         const E_Int posroe, const E_Float gamma,
                         FldArrayF& velo, FldArrayF& press, FldArrayF& mach)
{
  if (mach.getSize() != 0) return; // deja calcule
  
  E_Int npts = f.getSize(); mach.malloc(npts);
  computePressure(f, posro, posrou, posrov, posrow, posroe,
                  gamma, velo, press);
  E_Float u2, ainv;
  E_Float gammai = 1./gamma;
  E_Float* vx = velo.begin(1);
  E_Float* vy = velo.begin(2);
  E_Float* vz = velo.begin(3);
  const E_Float* rop = f.begin(posro);
  E_Float* pressp = press.begin();
  E_Float* machp = mach.begin();
  for (E_Int ind = 0; ind < npts; ind++)
  {
    u2 = vx[ind]*vx[ind] + vy[ind]*vy[ind] + vz[ind]*vz[ind];   
    ainv = gammai * rop[ind] / pressp[ind];
    machp[ind] = sqrt(u2 * ainv);
  }
}
