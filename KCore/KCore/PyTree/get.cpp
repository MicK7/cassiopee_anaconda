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
#include "PyTree/PyTree.h"
#include "String/kstring.h"
using namespace K_FLD;
using namespace std;

//==============================================================================
// Recherche par nom d'un seul niveau, retourne un seul noeud.
// IN: o: objet representant un noeud de pyTree
// IN: name: le nom du noeud
// OUT: retourne le nouveau noeud si trouve, sinon retourne NULL
//==============================================================================
PyObject* K_PYTREE::getNodeFromName1(PyObject* o, const char* name)
{
  PyObject* childrens = PyList_GetItem(o, 2);
  E_Int n = PyList_Size(childrens);
  PyObject *l; PyObject *node; char* str;
  for (E_Int i = 0; i < n; i++)
  {
    l = PyList_GetItem(childrens, i);
    node = PyList_GetItem(l, 0);
    str = PyString_AsString(node);
    if (K_STRING::cmp(str, name) == 0) return l;
  }
  return NULL;
}

//==============================================================================
// Recherche par nom d'un seul niveau, retourne une liste de noeuds.
// IN: o: objet representant un noeud de pyTree
// IN: type: le type du noeud
// OUT: out: la liste des noeuds trouves
//==============================================================================
void K_PYTREE::getNodesFromType1(PyObject* o, const char* type,
                                 vector<PyObject*>& out)
{
  PyObject* childrens = PyList_GetItem(o, 2);
  E_Int n = PyList_Size(childrens);
  PyObject *l; PyObject *node; char* str;
  for (E_Int i = 0; i < n; i++)
  {
    l = PyList_GetItem(childrens, i);
    node = PyList_GetItem(l, 3);
    str = PyString_AsString(node);
    if (K_STRING::cmp(str, type) == 0) out.push_back(l);
  }
}

//==============================================================================
// Retourne le nom du noeud
// Le tableau de char retourne est partage avec python.
// IN: o: noeud du pyTree
// OUT: retourne le ptr sur la chaine partagee
//==============================================================================
char* K_PYTREE::getNodeName(PyObject* o, vector<PyArrayObject*>& hook)
{
  PyObject* v = PyList_GetItem(o, 0);
  if (PyString_Check(v) == true)
  { char* r = PyString_AsString(v); return r; }
  return NULL;
}

//==============================================================================
// Retourne le type du noeud
// Le tableau de char retourne est partage avec python.
// IN: o: noeud du pyTree
// OUT: retourne le ptr sur la chaine partagee
//==============================================================================
char* K_PYTREE::getNodeType(PyObject* o, vector<PyArrayObject*>& hook)
{
  PyObject* v = PyList_GetItem(o, 3);
  if (PyString_Check(v) == true)
  { char* r = PyString_AsString(v); return r; }
  return NULL;
}

//==============================================================================
// Recupere la chaine de char si le noeud o contient une string ou un numpy char
// Le tableau de char retourne est partage avec python.
// Attention: dans le cas de numpy, il n'y a pas de \0 a la fin.
// IN: o: noeud du pyTree
// OUT: hook: contient les numpy char
// OUT: retourne le ptr sur la chaine partagee
//==============================================================================
char* K_PYTREE::getValueS(PyObject* o, vector<PyArrayObject*>& hook)
{
  //IMPORTNUMPY;
  PyObject* v = PyList_GetItem(o, 1);
  if (PyString_Check(v) == true) return PyString_AsString(v);
  PyArrayObject* ac = (PyArrayObject*)v; Py_INCREF(ac);
  char* d = (char*)PyArray_DATA(ac);
  hook.push_back(ac);
  return d;
}
//==============================================================================
// Recupere la chaine de char si le noeud o contient une string ou un numpy char
// Le tableau de char retourne est partage avec python.
// Attention: dans le cas de numpy, il n'y a pas de \0 a la fin.
// IN: o: noeud du pyTree
// OUT: s: taille de la chaine retournee
// OUT: hook: contient les numpy char
// OUT: retourne le ptr sur la chaine partagee
//==============================================================================
char* K_PYTREE::getValueS(PyObject* o, E_Int& s, vector<PyArrayObject*>& hook)
{
  IMPORTNUMPY;
  PyObject* v = PyList_GetItem(o, 1);
  if (PyString_Check(v) == true)
  { char* r = PyString_AsString(v); s = strlen(r); return r; }
  PyArrayObject* ac = (PyArrayObject*)v; Py_INCREF(ac);
  char* d = (char*)PyArray_DATA(ac);
  s = PyArray_DIM(ac, 0);
  hook.push_back(ac);
  return d;
}

//==============================================================================
// Retourne un pointeur sur un tableau d'entiers si le noeud o contient
// un numpy d'entiers
// IN: o: noeud de pyTree
// OUT: hook: contient les numpy d'entiers
// OUT: retourne le ptr sur le tableau d'entiers
//==============================================================================
E_Int* K_PYTREE::getValueAI(PyObject* o, vector<PyArrayObject*>& hook)
{
  //IMPORTNUMPY;
  PyObject* v = PyList_GetItem(o, 1);
  PyArrayObject* ac = (PyArrayObject*)v; Py_INCREF(ac);
  E_Int* d = (E_Int*)PyArray_DATA(ac);
  hook.push_back(ac);
  return d;
}

//==============================================================================
// Retourne un entier si le noeud o contient un numpy d'entiers ou un entier
// IN: o: noeud de pyTree
// OUT: retourne l'entier
//==============================================================================
E_Int K_PYTREE::getValueI(PyObject* o)
{
  //IMPORTNUMPY;
  PyObject* v = PyList_GetItem(o, 1);
  if (PyLong_Check(v) == true) return PyLong_AsLong(v);
  PyArrayObject* ac = (PyArrayObject*)v;
  E_Int* d = (E_Int*)PyArray_DATA(ac);
  E_Int ret = d[0];
  Py_DECREF(ac);
  return ret;
}

//==============================================================================
// Retourne un pointeur sur un tableau d'entier si le noeud o contient
// un numpy d'entiers
// IN: o: noeud de pyTree
// OUT: s0, s1: taille du tableau d'entiers (s0,s1)
// OUT: hook: contient les numpy d'entiers
// OUT: retourne le ptr sur le tableau d'entiers
//==============================================================================
E_Int* K_PYTREE::getValueAI(PyObject* o, E_Int& s0, E_Int& s1,
                            vector<PyArrayObject*>& hook)
{
  //IMPORTNUMPY;
  PyObject* v = PyList_GetItem(o, 1);
  PyArrayObject* ac = (PyArrayObject*)v; Py_INCREF(ac);
  E_Int* d = (E_Int*)PyArray_DATA(ac);
  int ndim = PyArray_NDIM(ac);
  if (ndim == 1) { s0 = PyArray_DIM(ac,0); s1 = 1; }
  else { s0 = PyArray_DIM(ac,0); s1 = PyArray_DIM(ac,1); }
  hook.push_back(ac);
  return d;
}

//==============================================================================
// Retourne un pointeur sur un tableau de doubles si le noeud o contient
// un numpy de doubles
// IN: o: noeud de pyTree
// OUT: hook: contient les numpy de doubles
// OUT: retourne le ptr sur le tableau de doubles
//==============================================================================
E_Float* K_PYTREE::getValueAF(PyObject* o, vector<PyArrayObject*>& hook)
{
  //IMPORTNUMPY
  PyObject* v = PyList_GetItem(o, 1);
  PyArrayObject* ac = (PyArrayObject*)v; Py_INCREF(ac);
  E_Float* d = (E_Float*)PyArray_DATA(ac);
  hook.push_back(ac);
  return d;
}

//==============================================================================
// Retourne un float si le noeud o contient un numpy de floats ou un float
// IN: o: noeud de pyTree
// OUT: retourne le float
//==============================================================================
E_Float K_PYTREE::getValueF(PyObject* o)
{
  IMPORTNUMPY;
  PyObject* v = PyList_GetItem(o, 1);
  if (PyFloat_Check(v) == true) return PyFloat_AsDouble(v);
  PyArrayObject* ac = (PyArrayObject*)v;
  E_Float* d = (E_Float*)PyArray_DATA(ac);
  E_Float ret = d[0];
  Py_DECREF(ac);
  return ret;
}

//==============================================================================
// fill name from pt
//==============================================================================
E_Int fillName(char*& pt, char* name)
{
  E_Int c = 0;
  while (*pt != '\0' && *pt != '/' && c < 255)
  { name[c] = *pt; c++; pt++; }
  return c;
}

//==============================================================================
// Recherche par path, retourne un seul noeud.
// IN: o: objet representant un noeud de pyTree (depart)
// IN: path: le chemin
// OUT: retourne le nouveau noeud si trouve, sinon retourne NULL
//==============================================================================
PyObject* K_PYTREE::getNodeFromPath(PyObject* o, const char* path)
{
  PyObject *l; PyObject *node; char* str;
  char name[256];
  char* pt = (char*)path;
  PyObject* next = o;
  E_Int c = 0;
  E_Boolean found;

  while (1)
  {
    // fill name
    c = 0;
    while (c == 0 && *pt != '\0')
    { c = fillName(pt, name); if (*pt == '/') pt++; }
    name[c] = '\0';
    // printf("name=%s\n", name);

    // find name
    PyObject* childrens = PyList_GetItem(next, 2);
    E_Int n = PyList_Size(childrens);

    found = false;
    for (E_Int i = 0; i < n; i++)
    {
      l = PyList_GetItem(childrens, i);
      node = PyList_GetItem(l, 0);
      str = PyString_AsString(node);
      // printf("comparing %s %s %c\n", str, name, *pt);
      if (K_STRING::cmp(str, name) == 0) { next = l; found = true; break; }
    }
    if (found == false) return NULL;
    if (*pt == '\0') return next;
  }
  return next;
}
