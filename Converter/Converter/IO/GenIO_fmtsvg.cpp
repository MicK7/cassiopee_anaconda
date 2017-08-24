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

// Formated svg file support

# include <stdio.h>
# include <stdlib.h>
# include <vector>
# include <string.h>

# include "GenIO.h"
# include "Array/Array.h"
# include "Def/DefFunction.h"
# include "Connect/connect.h"
# include "CompGeom/compGeom.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* svgread 
   Read svg polyLines and Bezier curves as i-arrays. 
   On regle la densite de points. */
//=============================================================================
E_Int K_IO::GenIO::svgread(
  char* file, char*& varString, E_Float density,
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, 
  vector<char*>& zoneNames)
{
  FldArrayF bezierPts(4, 3); // pts de controle Bezier cubique
  FldArrayF tmp;

  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: svgread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");

  // Lecture de l'entete
  E_Int res;
  char* buf = new char[BUFSIZE];
  char* prevData = new char[BUFSIZE];
  char* keyword = new char[BUFSIZE];
  E_Float zlayer = 0.;
  E_Float stack[3];
  E_Float xp0, yp0, zp0, xp1, yp1, zp1, alpha;

  res = readGivenKeyword(ptrFile, "<?");
  if (res == 0)
  {
    printf("Warning: svgread: cannot find xml header.\n");
    fclose(ptrFile);
    return 1;
  }
  res = readGivenKeyword(ptrFile, ">");
  if (res == 0)
  {
    printf("Warning: svgread: file seems to be corrupted.\n");
    fclose(ptrFile);
    return 1;
  }
  
  // Lecture du fichier
  list<char*> knownKeywords;
  knownKeywords.push_back("STYLE");
  knownKeywords.push_back("ID");
  knownKeywords.push_back("D");
  knownKeywords.push_back("TRANSFORM");

  res = readWord(ptrFile, buf);
  while (res >= 0)
  {
    // Commentaire
    if (strcmp(buf, "<!") == 0)
    {
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }
    // Entete svg
    else if (strcmp(buf, "<svg") == 0)
    {
      //printf("Entete svg\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Entete sopodi
    else if (strcmp(buf, "<sopodi:namedview") == 0)
    {
      //printf("Entete sopodi\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Entete metadata
    else if (strcmp(buf, "<metadata") == 0)
    {
      //printf("Bloc metadata\n");
      res = readGivenKeyword(ptrFile, "/METADATA>");
      if (res == 0) res = -1;
    }
    
    // Defs
    else if (strcmp(buf, "<defs") == 0)
    {
      //printf("Bloc defs\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Graphic layer
    else if (strcmp(buf, "<g") == 0)
    {
      zlayer = zlayer + 1.;
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Text
    else if (strcmp(buf, "<text") == 0)
    {
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Path => Objet
    else if (strcmp(buf, "<path") == 0)
    {
      FldArrayF* an = NULL;
      E_Int nmax = 1000;
      E_Int npts = 0;
      res = readDataAndKeyword(ptrFile, buf, knownKeywords, prevData, keyword);
      if (res == 1) res = -1;
      E_Boolean found = false;
      while (res >= 0 && found == false) // OK
      {
        if (strcmp(keyword, "D") == 0)
        {
          //an = new FldArrayF(nmax, 3);
          res = readWord(ptrFile, buf);
          E_Int reread = true;
          E_Int nctrlpts = 0;
          while (res >= 0)
          {
            E_Int lastChar = buf[strlen(buf)-1];
            if (lastChar == 'M' || lastChar == 'm')
            { 
              if (an != NULL)
              {
                an->reAllocMat(npts, 3);
                structField.push_back(an);
                ni.push_back(npts);
                nj.push_back(1);
                nk.push_back(1);
                npts = 0;
              }
              an = new FldArrayF(nmax, 3);
              reread = readTwoCoordinates(ptrFile, stack);
            }
            else if (lastChar == 'L' || lastChar == 'l') // lignes
            { 
              xp0 = stack[0]; yp0 = stack[1]; zp0 = zlayer;
              reread = readTwoCoordinates(ptrFile, stack);
              xp1 = stack[0]; yp1 = stack[1]; zp1 = zlayer;
              E_Float len = sqrt((xp1-xp0)*(xp1-xp0)+(yp1-yp0)*(yp1-yp0)+(zp1-zp0)*(zp1-zp0));
              E_Int NptsLine = E_Int(len*density)+2;
              alpha = 1. / (NptsLine-1);
              if (npts+NptsLine > nmax) 
              { nmax = nmax+2*NptsLine; an->reAllocMat(nmax, 3); }
              FldArrayF& coord = *an;

              for (E_Int i = 0; i < NptsLine-1; i++)
              {
                coord(npts, 1) = xp0 + i*alpha*(xp1 - xp0);
                coord(npts, 2) = yp0 + i*alpha*(yp1 - yp0);
                coord(npts, 3) = zp0 + i*alpha*(zp1 - zp0);
                npts++;
              }
            }
            else if (lastChar == 'C' || lastChar == 'c') // courbes
            {
              FldArrayF& coord = *an;
              bezierPts(nctrlpts, 1) = stack[0];
              bezierPts(nctrlpts, 2) = stack[1];
              bezierPts(nctrlpts, 3) = zlayer;
              nctrlpts++;
              for (E_Int i = 0; i < 3; i++)
              {
                reread = readTwoCoordinates(ptrFile, stack);
                bezierPts(nctrlpts, 1) = stack[0];
                bezierPts(nctrlpts, 2) = stack[1];
                bezierPts(nctrlpts, 3) = zlayer;
                nctrlpts++;
              }
              
              K_COMPGEOM::regularBezier(bezierPts.getSize(), -1, density,
                                        bezierPts.begin(1), 
                                        bezierPts.begin(2), 
                                        bezierPts.begin(3), tmp);
              E_Int NptsCurve = tmp.getSize();
              if (npts+ NptsCurve > nmax)
              { nmax = nmax+2*NptsCurve; an->reAllocMat(nmax,3); }
              for (E_Int i = 0; i < NptsCurve-1; i++)
              {
                coord(npts, 1) = tmp(i, 1);
                coord(npts, 2) = tmp(i, 2);
                coord(npts, 3) = tmp(i, 3); 
                npts++;
              }
              nctrlpts = 0;
            }
            else if (lastChar == 'Z' || lastChar == 'z')
            {
              reread = false;
            }
            if (reread == true) res = readWord(ptrFile, buf);
            else
            {
              FldArrayF& coord = *an;
              coord(npts, 1) = stack[0];
              coord(npts, 2) = stack[1];
              coord(npts, 3) = zlayer;
              npts++;
              res = -1;
            }
          }

          found = true;
          res = readGivenKeyword(ptrFile, ">"); // fin du bloc path
          if (res == 0) res = -1;
        }
        else if (strcmp(keyword, "STYLE") == 0)
        {
          res = readGivenKeyword(ptrFile, "\"");
          res = readGivenKeyword(ptrFile, "\"");
          if (res == 0) res = -1;
        }
        else if (strcmp(keyword, "ID") == 0)
        {
          res = readGivenKeyword(ptrFile, "\"");
          res = readGivenKeyword(ptrFile, "\"");
          if (res == 0) res = -1;
        }
        if (found == false)
        {
          res = readDataAndKeyword(ptrFile, buf, 
                                   knownKeywords, prevData, keyword);
          if (res == 1) res = -1;
        }
      } // bloc path
      
      if (found == true)
      {
        an->reAllocMat(npts, 3);
        structField.push_back(an);
        ni.push_back(npts);
        nj.push_back(1);
        nk.push_back(1);
      }    
    }
    res = readWord(ptrFile, buf);
  }
  // Cree les noms des zones
  for (unsigned int i = 0; i < structField.size(); i++)
  {
    char* zoneName = new char [128];
    sprintf(zoneName, "Zone%d",i);
    zoneNames.push_back(zoneName);
  }

  delete [] buf;
  delete [] prevData;
  delete [] keyword;
  fclose(ptrFile);
  return 0;
}
//=============================================================================
/* svgread 
   Read svg polyLines and Bezier curves as i-arrays. 
   On regle le nombre de pts par courbe. */
//=============================================================================
E_Int K_IO::GenIO::svgread(
  char* file, char*& varString, E_Int NptsCurve, E_Int NptsLine, 
  vector<FldArrayF*>& structField,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType, 
  vector<char*>& zoneNames)
{
  FldArrayF tmp(NptsCurve, 3);
  FldArrayF bezierPts(4, 3); // pts de controle Bezier cubique

  /* File Opening */
  FILE* ptrFile;
  ptrFile = fopen(file, "r");

  if (ptrFile == NULL)
  {
    printf("Warning: svgread: cannot open file %s.\n", file);
    return 1;
  }

  varString = new char [8];
  strcpy(varString, "x,y,z");

  // Lecture de l'entete
  E_Int res;
  char* buf = new char[BUFSIZE];
  char* prevData = new char[BUFSIZE];
  char* keyword = new char[BUFSIZE];
  E_Float zlayer = 0.;
  E_Float stack[3];
  E_Float xp0, yp0, zp0, xp1, yp1, zp1, alpha;

  res = readGivenKeyword(ptrFile, "<?");
  if (res == 0)
  {
    printf("Warning: svgread: cannot find xml header.\n");
    fclose(ptrFile);
    return 1;
  }
  res = readGivenKeyword(ptrFile, ">");
  if (res == 0)
  {
    printf("Warning: svgread: file seems to be corrupted.\n");
    fclose(ptrFile);
    return 1;
  }
  
  // Lecture du fichier
  list<char*> knownKeywords;
  knownKeywords.push_back("STYLE");
  knownKeywords.push_back("ID");
  knownKeywords.push_back("D");
  knownKeywords.push_back("TRANSFORM");

  res = readWord(ptrFile, buf);
  while (res >= 0)
  {
    // Commentaire
    if (strcmp(buf, "<!") == 0)
    {
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }
    // Entete svg
    else if (strcmp(buf, "<svg") == 0)
    {
      //printf("Entete svg\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Entete sopodi
    else if (strcmp(buf, "<sopodi:namedview") == 0)
    {
      //printf("Entete sopodi\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Entete metadata
    else if (strcmp(buf, "<metadata") == 0)
    {
      //printf("Bloc metadata\n");
      res = readGivenKeyword(ptrFile, "/METADATA>");
      if (res == 0) res = -1;
    }
    
    // Defs
    else if (strcmp(buf, "<defs") == 0)
    {
      //printf("Bloc defs\n");
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Graphic layer
    else if (strcmp(buf, "<g") == 0)
    {
      zlayer = zlayer + 1.;
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Text
    else if (strcmp(buf, "<text") == 0)
    {
      res = readGivenKeyword(ptrFile, ">");
      if (res == 0) res = -1;
    }

    // Path => Objet
    else if (strcmp(buf, "<path") == 0)
    {
      FldArrayF* an = NULL;
      E_Int nmax = NptsLine+NptsCurve;
      E_Int npts = 0;
      res = readDataAndKeyword(ptrFile, buf, knownKeywords, prevData, keyword);
      if (res == 1) res = -1;
      E_Boolean found = false;
      while (res >= 0 && found == false) // OK
      {
        if (strcmp(keyword, "D") == 0)
        {
          //an = new FldArrayF(nmax, 3);
          res = readWord(ptrFile, buf);
          E_Int reread = true;
          E_Int nctrlpts = 0;
          while (res >= 0)
          {
            E_Int lastChar = buf[strlen(buf)-1];
            if (lastChar == 'M' or lastChar == 'm')
            { 
              if (an != NULL)
              {
                an->reAllocMat(npts, 3);
                structField.push_back(an);
                ni.push_back(npts);
                nj.push_back(1);
                nk.push_back(1);
                npts = 0;
              }
              an = new FldArrayF(nmax, 3);
              reread = readTwoCoordinates(ptrFile, stack);
            }
            else if (lastChar == 'L' || lastChar == 'l') // lignes
            { 
              FldArrayF& coord = *an;
              if (npts+NptsLine > nmax) 
              { nmax = nmax+2*NptsLine; an->reAllocMat(nmax, 3); }
              xp0 = stack[0];
              yp0 = stack[1];
              zp0 = zlayer;
              reread = readTwoCoordinates(ptrFile, stack);
              xp1 = stack[0];
              yp1 = stack[1];
              zp1 = zlayer;
              alpha = 1. / (NptsLine-1);
              for (E_Int i = 0; i < NptsLine-1; i++)
              {
                coord(npts, 1) = xp0 + i*alpha*(xp1 - xp0);
                coord(npts, 2) = yp0 + i*alpha*(yp1 - yp0);
                coord(npts, 3) = zp0 + i*alpha*(zp1 - zp0);
                npts++;
              }
            }
            else if (lastChar == 'C' || lastChar == 'c') // courbes
            {
              FldArrayF& coord = *an;
              bezierPts(nctrlpts, 1) = stack[0];
              bezierPts(nctrlpts, 2) = stack[1];
              bezierPts(nctrlpts, 3) = zlayer;
              nctrlpts++;
              for (E_Int i = 0; i < 3; i++)
              {
                reread = readTwoCoordinates(ptrFile, stack);
                bezierPts(nctrlpts, 1) = stack[0];
                bezierPts(nctrlpts, 2) = stack[1];
                bezierPts(nctrlpts, 3) = zlayer;
                nctrlpts++;
              }
              K_COMPGEOM::bezier(bezierPts.getSize(), NptsCurve, 
                                 bezierPts.begin(1), 
                                 bezierPts.begin(2), 
                                 bezierPts.begin(3), tmp);

              if (npts+NptsCurve > nmax)
              { nmax = nmax+2*NptsCurve; an->reAllocMat(nmax,3); }
              for (E_Int i = 0; i < NptsCurve-1; i++)
              {
                coord(npts, 1) = tmp(i, 1);
                coord(npts, 2) = tmp(i, 2);
                coord(npts, 3) = tmp(i, 3); 
                npts++;
              }
              nctrlpts = 0;
            }
            else if (lastChar == 'Z' || lastChar == 'z')
            {
              reread = false;
            }
            if (reread == true)
              res = readWord(ptrFile, buf);
            else
            {
              FldArrayF& coord = *an;
              coord(npts, 1) = stack[0];
              coord(npts, 2) = stack[1];
              coord(npts, 3) = zlayer;
              npts++;
              res = -1;
            }
          }
          found = true;
          res = readGivenKeyword(ptrFile, ">"); // fin du bloc path
          if (res == 0) res = -1;
        }
        else if (strcmp(keyword, "STYLE") == 0)
        {
          res = readGivenKeyword(ptrFile, "\"");
          res = readGivenKeyword(ptrFile, "\"");
          if (res == 0) res = -1;
        }
        else if (strcmp(keyword, "ID") == 0)
        {
          res = readGivenKeyword(ptrFile, "\"");
          res = readGivenKeyword(ptrFile, "\"");
          if (res == 0) res = -1;
        }
        if (found == false)
          res = readDataAndKeyword(ptrFile, buf, 
                                   knownKeywords, prevData, keyword);
        if (res == 1) res = -1;
      } // bloc path
      
      if (found == true)
      {
        an->reAllocMat(npts, 3);
        structField.push_back(an);
        ni.push_back(npts);
        nj.push_back(1);
        nk.push_back(1);
      }
    }
    res = readWord(ptrFile, buf);
  }
  delete [] buf;
  delete [] prevData;
  delete [] keyword;
  fclose(ptrFile);
  return 0;
}

//=============================================================================
/* Write arrays as svg polyLines.
   Others are discarded. */
//=============================================================================
E_Int K_IO::GenIO::svgwrite(
  char* file, char* dataFmt, char* varString,
  vector<E_Int>& ni, vector<E_Int>& nj, vector<E_Int>& nk,
  vector<FldArrayF*>& structField,
  vector<FldArrayF*>& unstructField,
  vector<FldArrayI*>& connect,
  vector<E_Int>& eltType,   
  vector<char*>& zoneNames)
{
  // All zones must have posx, posy, posz
  E_Int posx, posy, posz, ind;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    printf("Warning: svgwrite: zones do not have coordinates. Not written.\n");
    return 1;
  }
  posx++; posy++; posz++;

  // Build writing data format
  char format1[20], format2[20], format3[20], format4[20];
  char dataFmtl[40];
  strcpy(dataFmtl, dataFmt);
  int l = strlen(dataFmt); 
  if (dataFmt[l-1] == ' ') dataFmtl[l-1] = '\0';

  // Build format{i}
  sprintf(format1," %s,%s", dataFmtl, dataFmtl);
  sprintf(format2," L %s,%s", dataFmtl, dataFmtl);
  sprintf(format3," %s,%s ", dataFmtl, dataFmtl);
  sprintf(format4," L %s,%s ", dataFmtl, dataFmtl);

  // Ecriture de l'entete
  FILE* ptrFile = fopen(file, "w");
  if (ptrFile == NULL) 
  {
    printf("Warning: svgwrite: I can't open file %s.\n", file);
    return 1;
  }

  fprintf(ptrFile, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
  fprintf(ptrFile, "<!-- Created with Converter (Onera) -->\n");

  // Entete svg
  fprintf(ptrFile, "<svg\n   xmlns:svg=\"http://www.w3.org/2000/svg\"\n   xmlns=\"http://www.w3.org/2000/svg\"\n   version=\"1.0\"\n   width=\"10.\"\n   height=\"10.\"\n   id=\"svg2\">\n");
  fprintf(ptrFile, "<defs\n   id=\"defs4\" />\n");
  fprintf(ptrFile, "<g\n   id=\"layer1\">\n");

  E_Int nzones = structField.size();
  E_Int nzoneu = unstructField.size();
  E_Int nc = 0;

  // Structured -> open polylines for each constant lines
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    E_Int nil = ni[zone];
    E_Int njl = nj[zone];
    E_Int nkl = nk[zone];
    FldArrayF& f = *structField[zone];

    for (E_Int k = 0; k < nkl; k++)
      for (E_Int j = 0; j < njl; j++)
      {
        fprintf(ptrFile, "<path\n   d=\"M");
        for (E_Int i = 0; i < nil; i++)
        {
          ind = i + j*nil + k*nil*njl;
          if (i == 0)
            fprintf(ptrFile, format1, f(ind,posx), f(ind,posy));
          else
            fprintf(ptrFile, format2, f(ind,posx), f(ind,posy));
        }
        fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
      }

    if (njl != 1)
    {
      for (E_Int k = 0; k < nkl; k++)
        for (E_Int i = 0; i < nil; i++)
        {
          fprintf(ptrFile, "<path\n   d=\"M");
          for (E_Int j = 0; j < njl; j++)
          {
            ind = i + j*nil + k*nil*njl;
            if (j == 0)
              fprintf(ptrFile,format1, f(ind,posx), f(ind,posy));
            else
              fprintf(ptrFile,format2, f(ind,posx), f(ind,posy));
          }
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
        }
    }

    if (nkl != 1)
    {
      for (E_Int j = 0; j < njl; j++)
        for (E_Int i = 0; i < nil; i++)
        {
          fprintf(ptrFile, "<path\n   d=\"M");
          for (E_Int k = 0; k < nkl; k++)
          {
            ind = i + j*nil + k*nil*njl;
            if (k == 0)
              fprintf(ptrFile,format1, f(ind,posx), f(ind,posy));
            else
              fprintf(ptrFile,format2, f(ind,posx), f(ind,posy));
          }
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
        }
    }
  }

  // Unstructured -> closed polylines
  for (E_Int zone = 0; zone < nzoneu; zone++)
  {
    FldArrayF& f = *unstructField[zone];
    FldArrayI& c = *connect[zone];

    for (E_Int i = 0; i < c.getSize(); i++)
    {
    
      switch (eltType[zone])
      {
        case 1: // BAR
          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) - 1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,2) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
          break;

        case 2: // TRI
          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) - 1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,2) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,3) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,1) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
          break;
          
        case 3: // QUAD
          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) - 1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,2) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,3) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,4) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,1) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
          break;

        case 4: // TETRA
          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) - 1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,2) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,3) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,1) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,4) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,2) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,4) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,3) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
          break;

        case 6: // PENTA (PRISM)
          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) - 1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,2) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,3) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,1) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,4) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,5) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,2) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,4) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,6) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,3) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,6) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,5) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
          break;

        case 7: // HEXA
          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) - 1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,2) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,3) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,4) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,1) - 1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,1) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,5) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,6) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,2) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,4) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,8) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,7) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          ind = c(i,3) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,5) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,8) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;

          fprintf(ptrFile, "<path\n   d=\"M");
          ind = c(i,6) -1;
          fprintf(ptrFile,format3, f(ind,posx), f(ind,posy));
          ind = c(i,7) -1;
          fprintf(ptrFile,format4, f(ind,posx), f(ind,posy));
          fprintf(ptrFile, "\"\nstyle=\"fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"\n   id=\"path%d\" />\n", nc); nc++;
          break;

        default:
          printf("Error: unrecognised element type.\n");
      }
    }
  }

  fprintf(ptrFile, "   </g>\n</svg>\n");
  fclose(ptrFile);
  return 0;
}
