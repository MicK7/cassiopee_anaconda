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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "Polygon.h"
#include "Connect/IdTool.h"

#include <vector>
#define Vector_t std::vector

#ifdef DEBUG_POLYGON
#include "IO/DynArrayIO.h"
#include <iostream>
#endif

namespace K_MESH
{ 

//=============================================================================
/*void Polygon::getBoundary
(const Polygon&  T1, const Polygon&  T2, K_MESH::NO_Edge& b)
{
  K_MESH::NO_Edge Ei, Ej;
  for (E_Int i = 0; i < NB_NODES; ++i)
  {
    T1.getBoundary(i, Ei);
    for (E_Int j = 0; j < NB_NODES; ++j)
    {
      T2.getBoundary(j, Ej);
      if (Ei == Ej)
      {
        b = Ei; return;
      }
    }
  }
}
*/
//=============================================================================
void Polygon::getBoundary
(const Polygon&  T1, const Polygon&  T2, E_Int& i1, E_Int& i2)
{
  K_MESH::NO_Edge Ei, Ej;
  for (i1=0; i1 < T1.NB_NODES; ++i1)
  {
    T1.getBoundary(i1, Ei);
    for (i2=0; i2 < T2.NB_NODES; ++i2)
    {
      T2.getBoundary(i2, Ej);
      if (Ei == Ej)
        return;
    }
  }
}

//=============================================================================
E_Int Polygon::getOrientation
(const Polygon&  PG, const E_Int& Ni, const E_Int& Nj, E_Bool& same_orient)
{
  same_orient = false;
 
  for (E_Int n = 0; n < PG.NB_NODES; ++n)
  {
    if ((PG._nodes[n] == Ni) && (PG._nodes[(n+1)%PG.NB_NODES] == Nj))
    {
      same_orient = true;
      return 0;
    }
    if ((PG._nodes[n] == Nj) && (PG._nodes[(n+1)%PG.NB_NODES] == Ni))
      return 0;
  }

  return -1;
}

///
bool Polygon:: is_convex
(const K_FLD::FloatArray& coord, const E_Int* nodes, E_Int nb_nodes, E_Int index_start,
const E_Float* normal, E_Float convexity_tol, E_Int& iworst)
{

  E_Float worst_ps = K_CONST::E_MAX_FLOAT, Ei[3], Ej[3], ps;
  bool convex = true;
  iworst = E_IDX_NONE;

  E_Float Z[3];
  for (E_Int i = 1; i < nb_nodes + 1; ++i)
  {
    E_Int ei = nodes[i%nb_nodes] - index_start;
    E_Int eim1 = nodes[i - 1] - index_start;
    E_Int eip1 = nodes[(i + 1) % nb_nodes] - index_start;

    K_FUNC::diff<3>(coord.col(ei), coord.col(eim1), &Ei[0]); //fixme : no normalization ??
    K_FUNC::diff<3>(coord.col(eip1), coord.col(ei), &Ej[0]);
    K_FUNC::sum<3>(normal, coord.col(ei), Z);

    E_Float uu = K_FUNC::zzdet4(coord.col(eim1), coord.col(ei), coord.col(eip1), Z);

    if (uu < -convexity_tol) //concavity
    {
      convex = false;
      ps = K_FUNC::dot<3>(Ei, Ej);
      if (ps < worst_ps)
        iworst = i%nb_nodes;
    }
  }

  return convex;
}

///
bool Polygon::is_convex
(const K_FLD::FloatArray& coord, const E_Int* nodes, E_Int nb_nodes, E_Int index_start, const E_Float *normal, E_Float convexity_tol)
{
  E_Float Z[3], Ei[3], Ej[3];
  for (E_Int i = 1; i < nb_nodes + 1; ++i)
  {
    E_Int ei = nodes[i%nb_nodes] - index_start;
    E_Int eim1 = nodes[i - 1] - index_start;
    E_Int eip1 = nodes[(i + 1) % nb_nodes] - index_start;

    K_FUNC::diff<3>(coord.col(ei), coord.col(eim1), &Ei[0]);
    K_FUNC::diff<3>(coord.col(eip1), coord.col(ei), &Ej[0]);

    /*K_FUNC::normalize<3>(Ei);
    K_FUNC::normalize<3>(Ej);
    K_FUNC::crossProduct<3>(Ei, Ej, Ek);

    ps = K_FUNC::dot<3>(Ek, normals.col(globalId[K]));
    */
    K_FUNC::sum<3>(normal, coord.col(ei), Z);

    E_Float uu = K_FUNC::zzdet4(coord.col(eim1), coord.col(ei), coord.col(eip1), Z);
    /*
    if ((ps*uu < 0.) && (::fabs(ps) > E_EPSILON) && (::fabs(uu) > E_EPSILON))
    {
    assert(false);
    }*/

    if (uu < -convexity_tol) //concavity
      return false;
  }

  return true;
}


}
