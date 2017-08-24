/*    
    Copyright 2013-2016 Onera.

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
#ifndef __K_MESH_POLYHEDRON_H__
#define __K_MESH_POLYHEDRON_H__

#include "MeshElement/Polygon.h"
#include "Fld/DynArray.h"
#include "Fld/ArrayAccessor.h"
#include "Fld/ngon_t.hxx"
#include "MeshElement/Tetrahedron.h"
#include "MeshElement/Polyhedron.h"
#include "Connect/GeomAlgo.h"

namespace K_MESH
{
  
//enum TopoShape { UNKNOWN, STAR_SHAPED, CONVEX, CONCAVE};
#define UNKNOWN 0
#define STAR_SHAPED 7

template <int TopoShape>
class Metric
{
  public:
  typedef ngon_t<K_FLD::IntArray> ngon_type;
  
  template<typename TriangulatorType>
  static inline E_Int volume(const TriangulatorType& t, const K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, E_Float&V, E_Float* G);
  
  static inline E_Float volumeTH4(const E_Float* p0, const E_Float* p1, const E_Float* p2, const E_Float* p3){return Tetrahedron::volume(p0,p1,p2,p3);}
  
  template<typename TriangulatorType>
  static inline E_Int volumes(const K_FLD::FloatArray& crd, const ngon_type& ng, std::vector<E_Float>& vols);
  
  template<typename TriangulatorType>
  static inline E_Int centroids(const ngon_type& ng, const K_FLD::FloatArray& crd, K_FLD::FloatArray& centroids);
  
  static inline void reorient_PHT3(K_FLD::IntArray& connectT3)
  {
    //Dual graph
    K_FLD::IntArray neighbors;
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectT3, neighbors);
  
    // Make facets orientation consistent
    std::vector<E_Int> orient;
    K_CONNECT::EltAlgo<K_MESH::Triangle>::reversi_connex(connectT3, neighbors, 0/*Kseed*/, orient);
    //permut the third and second nodes
    for (size_t k = 0; k < orient.size(); ++k)
      if (orient[k] == -1)
        std::swap(connectT3(1,k), connectT3(2,k));
    return;
  }
  
  static inline void reorient_PHT3(K_FLD::IntArray& connectT3,  K_FLD::IntArray& neighbors)
  {
    //Dual graph
    neighbors.clear();
    K_CONNECT::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(connectT3, neighbors);
  
    // Make facets orientation consistent
    std::vector<E_Int> orient;
    K_CONNECT::EltAlgo<K_MESH::Triangle>::reversi_connex(connectT3, neighbors, 0/*Kseed*/, orient);
    //permut the third and second nodes
    for (size_t k = 0; k < orient.size(); ++k)
      if (orient[k] == -1)
      {
        std::swap(connectT3(1,k), connectT3(2,k));
        std::swap(neighbors(1,k), neighbors(2,k));
      }
    return;
  }
  
  
  
  ///
  static E_Int is_concave
  (const K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, const ngon_unit& orient, bool & concave, E_Float threshold = E_EPSILON)
  {
    concave = false;
    
    // angular criterion
    E_Float angle_threshold = K_CONST::E_PI*(1. - threshold);
    angle_threshold = std::min(K_CONST::E_PI, angle_threshold);
    angle_threshold = std::max(angle_threshold, E_EPSILON);
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
  
    const E_Int & nb_polys = ng.PHs.stride(PHi);
  
    if (nb_polys == 4) // TH4
      return 0;
    
    // build PG neighborhood for that given PH
    Vector_t<E_Int> oids;
    ngon_unit lneighbors, lpgs;
    ng.PGs.extract(ng.PHs.get_facets_ptr(PHi), nb_polys, lpgs, oids);
    
    assert(lpgs.size() == nb_polys);
    
    ngon_type::build_pg_neighborhood(lpgs, lneighbors);
    
    E_Float ni[3], nj[3];

    for (E_Int PGi=0; PGi <nb_polys; ++PGi)
    {
      //std::cout << PGi << std::endl;
      E_Int* pNi      = lpgs.get_facets_ptr(PGi);
      E_Int  nb_nodes = lpgs.stride(PGi);
      E_Int* pKn      = lneighbors.get_facets_ptr(PGi);
      
      K_MESH::Polygon::normal<acrd_t, 3>(acrd, pNi, nb_nodes, 1, ni);
      bool reversed=false;
      if (orient.get_facet(PHi, PGi) == -1)
      {
        ni[0]=-ni[0];
        ni[1]=-ni[1];
        ni[2]=-ni[2];
        reversed=true;
      }
      
      for (E_Int j=0; j < nb_nodes; ++j)
      {
        E_Int e0 = *(pNi+j)-1;
        E_Int e1 = *(pNi+(j+1)%nb_nodes)-1;
        
        if (reversed)
          std::swap(e0,e1);
        
        const E_Float* E0 = crd.col(e0);
        const E_Float* E1 = crd.col(e1);
        
        E_Int PGj=*(pKn+j);
        if (PGj == E_IDX_NONE)
          return 1;
        
        const E_Int* pNj = lpgs.get_facets_ptr(PGj);
        E_Int nb_nodsj   = lpgs.stride(PGj);
        
        K_MESH::Polygon::normal<acrd_t, 3>(acrd, pNj, nb_nodsj, 1, nj);

        assert (PGj >= 0 && PGj < nb_polys);
         
        if (orient.get_facet(PHi, PGj) == -1)
        {
          nj[0]=-nj[0];
          nj[1]=-nj[1];
          nj[2]=-nj[2];
        }
                
        E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Polygon>::angle_measure(ni, nj, E0, E1);
        
        if (alpha < angle_threshold)
        {
          concave=true;
#ifdef DEBUG_POLYHEDRON
          K_FLD::FloatArray crdt(crd);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, lpgs, PGi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, lpgs, PGj);
          
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG1o.mesh", crd, lpgs, PGi, ni);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG2o.mesh", crd, lpgs, PGj, nj);
#endif
          return 0;
        }
      }
    }
  
    return 0;
  }
  
  ///
  static E_Int is_concave
  (const K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi, const ngon_unit& orient, bool & concave, 
   std::set<K_MESH::NO_Edge>& relex_edges, E_Float threshold = E_EPSILON)
  {
    concave = false;
    
    // angular criterion
    E_Float angle_threshold = K_CONST::E_PI*(1. - threshold);
    angle_threshold = std::min(K_CONST::E_PI, angle_threshold);
    angle_threshold = std::max(angle_threshold, E_EPSILON);
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t acrd(crd);
  
    const E_Int & nb_polys = ng.PHs.stride(PHi);
  
    if (nb_polys == 4) // TH4
      return 0;
    
    // build PG neighborhood for that given PH
    Vector_t<E_Int> oids;
    ngon_unit lneighbors, lpgs;
    ng.PGs.extract(ng.PHs.get_facets_ptr(PHi), nb_polys, lpgs, oids);
    
    assert(lpgs.size() == nb_polys);
    
    ngon_type::build_pg_neighborhood(lpgs, lneighbors);
    
    E_Float ni[3], nj[3];

    for (E_Int PGi=0; PGi <nb_polys; ++PGi)
    {
      //std::cout << PGi << std::endl;
      E_Int* pNi      = lpgs.get_facets_ptr(PGi);
      E_Int  nb_nodes = lpgs.stride(PGi);
      E_Int* pKn      = lneighbors.get_facets_ptr(PGi);
      
      K_MESH::Polygon::normal<acrd_t, 3>(acrd, pNi, nb_nodes, 1, ni);
      bool reversed=false;
      if (orient.get_facet(PHi, PGi) == -1)
      {
        ni[0]=-ni[0];
        ni[1]=-ni[1];
        ni[2]=-ni[2];
        reversed=true;
      }
      
      for (E_Int j=0; j < nb_nodes; ++j)
      {
        E_Int e0 = *(pNi+j)-1;
        E_Int e1 = *(pNi+(j+1)%nb_nodes)-1;
        
        if (reversed)
          std::swap(e0,e1);
        
        const E_Float* E0 = crd.col(e0);
        const E_Float* E1 = crd.col(e1);
        
        E_Int PGj=*(pKn+j);
        if (PGj == E_IDX_NONE)
          return 1;
        
        const E_Int* pNj = lpgs.get_facets_ptr(PGj);
        E_Int nb_nodsj   = lpgs.stride(PGj);
        
        K_MESH::Polygon::normal<acrd_t, 3>(acrd, pNj, nb_nodsj, 1, nj);

        assert (PGj >= 0 && PGj < nb_polys);
         
        if (orient.get_facet(PHi, PGj) == -1)
        {
          nj[0]=-nj[0];
          nj[1]=-nj[1];
          nj[2]=-nj[2];
        }
                
        E_Float alpha = K_CONNECT::GeomAlgo<K_MESH::Polygon>::angle_measure(ni, nj, E0, E1);
        
        if (alpha < angle_threshold)
        {
          concave=true;
          relex_edges.insert(K_MESH::NO_Edge(e0,e1));
#ifdef DEBUG_POLYHEDRON
          K_FLD::FloatArray crdt(crd);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, lpgs, PGi);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PG(crdt, lpgs, PGj);
          
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG1o.mesh", crd, lpgs, PGi, ni);
          NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_wired_PG("PG2o.mesh", crd, lpgs, PGj, nj);
#endif
        }
      }
    }
  
    return 0;
  }
  
};

///
template <int TopoShape>
class Polyhedron : public Metric<TopoShape>
{  
  
public:
  //static const E_Int NB_NODES;
  //static const E_Int NB_TRIS;

  typedef K_MESH::Polygon boundary_type;
    
public:
  Polyhedron(){}
  
  
  ///
  template <typename CoordAcc>
  static inline void iso_barycenter(const CoordAcc& coord, const E_Int* nodes, E_Int nb_nodes, E_Int index_start, E_Float* G)
  {
    //
    for (size_t d = 0; d < 3; ++d) G[d] = 0.;

    for (E_Int i = 0; i < nb_nodes; ++i)
    {
      for (size_t d = 0; d < 3; ++d)
      {
        //std::cout << "v : " << coord.getVal(node(i), d) << std::endl;
        G[d] += coord.getVal(nodes[i]-index_start, d);
      }
    }

    E_Float k = 1. / (E_Float)nb_nodes;

    for (size_t i = 0; i < 3; ++i) G[i] *= k;
    //std::cout << "G : " << G[0] << "/" << G[1] << "/" << G[2] << std::endl;
  }
  
private:
  
  Polyhedron(const Polyhedron& orig);

};

//SPECIALIZATIONS //////////////////////////////////////////////////////
/// not necessary to reorient since the centroid is inside so all the volume are consistent : either all positive or all negative
template<>
inline void K_MESH::Metric<STAR_SHAPED>::reorient_PHT3(K_FLD::IntArray& connectT3){return;}

//must be absolute value when not reorienting as for the star-shaped case
template<>
inline E_Float K_MESH::Metric<STAR_SHAPED>::volumeTH4(const E_Float* p0, const E_Float* p1, const E_Float* p2, const E_Float* p3)
{return ::fabs(Tetrahedron::volume(p0,p1,p2,p3));}
///////////////////////////////////////////////////////////////////////////

///
template<int TopoShape>
template<typename TriangulatorType>
E_Int K_MESH::Metric<TopoShape>::volume
(const TriangulatorType& t, const K_FLD::FloatArray& crd, const ngon_type& ng, E_Int PHi , E_Float&V, E_Float* G)
{
  G[0]=G[1]=G[2]=0.;
  
  V=K_CONST::E_MAX_FLOAT;
  E_Float p4[3], v;
  typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
  acrd_t acrd(crd);
  
  std::vector<E_Int> phnodes;
  bool is_th4=false;
  
  const E_Int & nb_faces = ng.PHs.stride(PHi);
  
  if (nb_faces == 4) // TH4
  {
    ng.unique_nodes_ph(PHi, phnodes);
    const E_Int& nb_nodes = phnodes.size();
    if (nb_nodes == 4)
      is_th4=true;
  }
  
  if (is_th4)
  {
    const E_Float* p1 = crd.col(phnodes[0]-1);
    const E_Float* p2 = crd.col(phnodes[1]-1);
    const E_Float* p3 = crd.col(phnodes[2]-1);
    const E_Float* p4 = crd.col(phnodes[3]-1);
    
    V=::fabs(::K_MESH::Tetrahedron::volume(p1, p2, p3, p4));
    
    G[0] = 0.25 * (p1[0]+p2[0]+p3[0]+p4[0]);
    G[1] = 0.25 * (p1[1]+p2[1]+p3[1]+p4[1]);
    G[2] = 0.25 * (p1[2]+p2[2]+p3[2]+p4[2]);
    
    return 0;
  }
  else if (phnodes.empty())
    ng.unique_nodes_ph(PHi, phnodes);//ng.nodes_ph(PHi, phnodes);
  
  K_MESH::Polyhedron<TopoShape>::iso_barycenter(acrd, &phnodes[0], phnodes.size(), 1, p4);
  
  K_FLD::IntArray connectT3;
  E_Int err = ngon_type::triangulate_ph(t, ng, PHi, crd, connectT3); // PH -> PHT3
  if (connectT3.cols() == 0 || err)
  {
    //std::cout << "could not triangulate properly" << std::endl;
    return 1;
  }
  
  //MIO::write("PHT3.mesh", crd, connectT3, "TRI");

  reorient_PHT3(connectT3); // does nothing if star shaped
  
  V=0.;
  for (size_t i=0; i < connectT3.cols(); ++i)
  {
    const K_FLD::IntArray::const_iterator pS = connectT3.col(i);
    const E_Float* p1 = crd.col(*pS);
    const E_Float* p2 = crd.col(*(pS+1));
    const E_Float* p3 = crd.col(*(pS+2));
    
    v=volumeTH4(p1, p2, p3, p4); //must be absolute value when not reorienting as for the star-shaped case
    
    V += v;
    
    G[0] += 0.25 * v * (p1[0]+p2[0]+p3[0]+p4[0]);
    G[1] += 0.25 * v * (p1[1]+p2[1]+p3[1]+p4[1]);
    G[2] += 0.25 * v * (p1[2]+p2[2]+p3[2]+p4[2]);
    
  }
  
  G[0] /= V;
  G[1] /= V;
  G[2] /= V;

  return 0;
}


///
template<int TopoShape>
template<typename TriangulatorType>
E_Int K_MESH::Metric<TopoShape>::volumes
(const K_FLD::FloatArray& crd, const ngon_type& ng, std::vector<E_Float>& vols)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();
  E_Int nb_phs = ng.PHs.size();
  //std::cout << "nb of pHs : " << nb_phs << std::endl;
  
  vols.clear();
  vols.resize(nb_phs, 0.);
  
  E_Float Gdum[3];
  TriangulatorType dt;
  
  for (size_t i=0; i < nb_phs; ++i)
    volume<TriangulatorType>(dt, crd, ng, i, vols[i], Gdum);
  
  return 0;
}

///
template<int TopoShape>
template<typename TriangulatorType>
E_Int K_MESH::Metric<TopoShape>::centroids(const ngon_type& ng, const K_FLD::FloatArray& crd, K_FLD::FloatArray& centroids)
{
  ng.PGs.updateFacets();
  ng.PHs.updateFacets();
  
  E_Int nb_phs = ng.PHs.size();
  //std::cout << "nb of pHs : " << nb_phs << std::endl;
  
  centroids.clear();
  centroids.resize(3, nb_phs, 0.);
  
  E_Float v;
  TriangulatorType dt;
  
  for (size_t i=0; i < nb_phs; ++i)
    volume<TriangulatorType>(dt, crd, ng, i, v, centroids.col(i));
  
  return 0;
}

}
#endif	/* __K_MESH_HEXAHEDRON_H__ */
