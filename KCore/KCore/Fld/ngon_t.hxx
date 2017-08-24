/*    
    Copyright 2013-2015 Onera.

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

#ifndef __NGON_T_HXX__
#define	__NGON_T_HXX__

#include "ngon_unit.h"

#include <iostream>
#include <map>
#include "MeshElement/Edge.h"
#include "MeshElement/Polygon.h"
#include "MeshElement/Triangle.h"
#include "MeshElement/Quadrangle.h"
#include "Connect/IdTool.h"
#include "Connect/EltAlgo.h"
#include "Connect/BARSplitter.h"

#ifdef DEBUG_NGON_T
#include "IO/io.h"
#include "Nuga/Boolean/NGON_debug.h"
#include "Nuga/Boolean/TRI_debug.h"
#define NGDBG NGON_debug<K_FLD::FloatArray, Connectivity_t>
#endif
#ifdef FLAG_STEP
#include "chrono.h"
#endif

#define STARTFACE 0
#define STARTELTS(cn) (2+cn[1])

typedef E_Int (*pg_transfo_func)(const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& tri_pgs, Vector_t<E_Int>& oids, E_Float optional_tol, const Vector_t<bool>* to_process);

template <typename Connectivity_t>
struct ngon_t
{
  ///
  ngon_t(const Connectivity_t& cNGON) :PGs(cNGON.begin() + STARTFACE), PHs(cNGON.begin() + STARTELTS(cNGON)){ PGs.updateFacets(); PHs.updateFacets(); }
  ///
  ngon_t(const ngon_t& ngt) :PGs(ngt.PGs), PHs(ngt.PHs){ PGs.updateFacets(); PHs.updateFacets(); }
  ///
  ngon_t& operator=(const Connectivity_t& cNGON)
  {
    PGs = ngon_unit(cNGON.begin() + STARTFACE); PGs.updateFacets();
    PHs = ngon_unit(cNGON.begin() + STARTELTS(cNGON)); PHs.updateFacets();
    return *this;
  }
  ///
  ngon_t(const ngon_unit& pg, const ngon_unit& ph) :PGs(pg), PHs(ph){ PGs.updateFacets(); PHs.updateFacets(); }
  ///
  ngon_t(const ngon_unit& pgs, bool one_ph_for_all=false):PGs(pgs)//false means one ph per pg.
  {
    PGs.updateFacets();
    E_Int sz= PGs.size();
        
    if (one_ph_for_all)
    {
      PHs._NGON.reserve(sz+3);
      PHs._NGON.resize(2,0);
      PHs._NGON.push_back(sz);
      for (E_Int i = 0; i < sz ; ++i) PHs._NGON.push_back(i+1);
      PHs._NGON[0] = 1;
      PHs._type.resize(1, INITIAL_SKIN);//only one PH => skin type
    }
    else
    {
      PHs._NGON.resize(2*sz+2, 1);
      for (E_Int i=0; i <sz ; ++i) PHs._NGON[2+2*i+1]=i+1;
      PHs._NGON[0] = sz;
    }
         
    PHs._NGON[1] = PHs._NGON.size()-2;
    PHs.updateFacets();
  }
  
  ///
  ngon_t(){}

  bool is_consistent(){ return PHs.is_consistent() && PGs.is_consistent(); }
  
  ///
  void addPH(const ngon_t& ngi, E_Int PHi)
  {
    if (PGs.size()) PGs.updateFacets();
    if (PHs.size()) PHs.updateFacets();
    
    E_Int npgid = PGs.size()+1;
    E_Int nb_pgs = ngi.PHs.stride(PHi);
    
    Vector_t<E_Int> molecPH;
        
    const E_Int* pPGi = ngi.PHs.get_facets_ptr(PHi);
    for (E_Int i = 0; i <nb_pgs; ++i, ++npgid)
    {
      E_Int PGi = *(pPGi+i)-1;
      PGs.add(ngi.PGs.stride(PGi), ngi.PGs.get_facets_ptr(PGi));
      if (!ngi.PGs._type.empty())
        PGs._type.push_back(ngi.PGs._type[PGi]);
      if (ngi.PGs._ancEs.cols() != 0)
        PGs._ancEs.pushBack(ngi.PGs._ancEs.col(PGi), ngi.PGs._ancEs.col(PGi)+2);

      molecPH.push_back(npgid);
    }
    
    PHs.add(molecPH.size(), &molecPH[0]);
    if (!ngi.PHs._type.empty())
      PHs._type.push_back(ngi.PHs._type[PHi]);
    if (ngi.PHs._ancEs.cols() != 0)
      PHs._ancEs.pushBack(ngi.PHs._ancEs.col(PHi), ngi.PHs._ancEs.col(PHi) + 2);
  }
  
  ///
  void export_to_array (Connectivity_t& c) const { 
    c.clear();
    //if (!PGs.size())
    //  return;
    //
    //c.pushBack(PGs._NGON);
    if (PGs.size())
      c.pushBack(PGs._NGON);
    else
      c.resize(1, 2, 0);

    if (PHs.size())
      c.pushBack(PHs._NGON);
    else
    {
      const E_Int& sizeFN = c[1];
      c.resize(1, sizeFN+4, 0);
    }
  }
  
  // Converts a SURFACE formatted as Edges/PGs into a ngon_unit (PHs is cleared) representing the PGs
  // WARNING : inconsistent orientation of the surface upon exit
  void export_surfacic_view(Connectivity_t& c)
  {
    c.clear();
    
    PGs.updateFacets();
    PHs.updateFacets();

    // Surfacic ngon ?
    for (E_Int i=0; (i<PGs.size()); ++i){
      if (PGs.stride(i) != 2)
        return;
    }
    
    c.resize(1, 2, 0);
     
    E_Int nb_polys=0;
    std::vector<E_Int> molec,snodes;
    K_FLD::IntArray cB;
        
    ngon_unit& polys = PHs;
    ngon_unit& edges = PGs;
    
    for (E_Int i=0; (i<polys.size()); ++i){
            
      E_Int nb_edges = polys.stride(i);
      const E_Int* pEs = polys.get_facets_ptr(i);
      
      cB.clear();
      for (size_t j=0; j < nb_edges; ++j)
      {
        E_Int Ej = pEs[j]-1;
        const E_Int* pN = edges.get_facets_ptr(Ej);
        cB.pushBack(pN, pN+2);
      }
      
      molec.clear();
      molec.push_back(nb_edges);// stride : nb_edges = nb nodes
      
      // sort the nodes
      BARSplitter::getSortedNodes(cB, snodes);
      
      /*if (snodes.size() != nb_edges) //degen
      {
        std::cout << cB << std::endl;
        continue;
      }*/
      
      ++nb_polys;
      
      molec.insert(molec.end(), snodes.begin(), snodes.end());
      
      //std::cout << molec.size() << std::endl;
      
      c.pushBack(molec);
    }
    
    c[0]=nb_polys;
    c[1]=c.getSize()-2;
    
  }
  
  // Converts a SURFACE formatted as Edges/PGs into a ngon_unit (PHs is cleared) representing the PGs
  // WARNING : inconsistent orientation of the surface upon exit
  // Morse version
  void export_surfacic_view(Connectivity_t& FN, Connectivity_t& xFN)
  {
    std::vector<E_Int> fn, xfn;

    const ngon_unit& polys = PHs;
    const ngon_unit& edges = PGs;
    
    edges.updateFacets();
    polys.updateFacets();
    
    E_Int nb_edges = edges.size();
    E_Int nb_pgs   = polys.size();

    // Surfacic ngon ?
    for (E_Int i=0; (i<nb_edges); ++i){
      if (edges.stride(i) != 2)
        return;
    }
    
    xfn.resize(nb_pgs+1);
    xfn[0]=0;
    
    //
    std::vector<E_Int> snodes;
    K_FLD::IntArray cB;

    //
    for (E_Int i=0; i<nb_pgs; ++i){
            
      E_Int nb_es = polys.stride(i);
      const E_Int* pEs = polys.get_facets_ptr(i);
      
      cB.clear();
      for (size_t j=0; j < nb_es; ++j)
      {
        E_Int Ej = pEs[j]-1;
        const E_Int* pN = edges.get_facets_ptr(Ej);
        cB.pushBack(pN, pN+2);
      }
                  
      // sort the nodes
      BARSplitter::getSortedNodes(cB, snodes);
      
      //if (snodes.size() != nb_es) //degen
        //continue;

      fn.insert(fn.end(), snodes.begin(), snodes.end());
      xfn[i+1]=fn.size();
    }
    
    FN = Connectivity_t(fn.size(), 1, &fn[0]);
    xFN = Connectivity_t(xfn.size(), 1, &xfn[0]);
  }

  ///
  inline bool empty() const {return (PGs._NGON.empty() || PHs._NGON.empty());}
  ///
  inline virtual void clear() {PGs.clear(); PHs.clear();}
  ///
  inline const E_Int* begin() const {return PGs.begin();}
    
  ///
  static void select_phs(const ngon_t& NGin, const Vector_t<bool>& keep, Vector_t<E_Int>& new_pg_ids, ngon_t& NGout)
  {
    new_pg_ids.clear();
    
    Vector_t<E_Int> stacked_pg_ids;
    __get_selected_ids(NGin.PHs, keep, new_pg_ids, stacked_pg_ids);
    __create_selected_phs(NGin.PHs, keep, new_pg_ids, NGout.PHs);
    __create_selected_pgs(NGin.PGs, stacked_pg_ids, NGout.PGs);
  }
  
  ///
  static void split_phs(const ngon_t& NGin, const Vector_t<bool>& flag, ngon_t& cng_ok, ngon_t& cng_nok)
  {
    ngon_t tmp(NGin); //use tmp because cng_ok or cng_nok could be NGin
    Vector_t<bool> nflag(flag); //copy at the beginnning because flag might change after compress (if it is NGin._external)
    K_CONNECT::IdTool::negative(nflag);
    Vector_t<E_Int> new_pg_ids;
    select_phs(tmp, flag, new_pg_ids, cng_ok);
    select_phs(tmp, nflag, new_pg_ids, cng_nok);
  }
  
  ///
  virtual void append(const ngon_t& NGin)
  {     
    if (NGin.empty())
      return;
    
    E_Int shiftv = PGs.size(); //to apply to the PHin
    
    // PGs    
    PGs.append(NGin.PGs);
    // PHs
    ngon_unit tmp(NGin.PHs);
    tmp.shift(shiftv);
    PHs.append(tmp);
  }
  
  ///
  E_Int flag_externals(E_Int FLAG)
  {
    if (PHs.size() == PGs.size()) // if one pg per ph : it's a surface so all PGs/PHs are skin ones.
    {
      PHs._type.resize(PHs.size(), FLAG);
      PGs._type.resize(PGs.size(), FLAG);
    }
    else
    {
      flag_external_pgs(FLAG);
      flag_external_phs(FLAG);
    }
  
    return 0;
  }
  
  ///
  E_Int flag_external_pgs(E_Int FLAG)
  {
    PHs.updateFacets();
    
    E_Int nb_faces = PGs.size();
    E_Int nb_elts = PHs.size();
    Vector_t<E_Int> face_count(nb_faces, 0);
  
    // Loop through the elements and increment face_count
    for (E_Int i = 0; i < nb_elts; ++i)
    {
      E_Int nb_facets = PHs.stride(i);
      for (E_Int j = 0; j < nb_facets; ++j)
      {
        const E_Int& Fi = PHs.get_facet(i,j);
        //std::cout << "Fi : " << Fi << std::endl;
         ++face_count[Fi-1];
      }
    }
  
    // External faces are those with a count equal to 1.
    PGs._type.resize(nb_faces, INNER);
    for (E_Int i = 0; i < nb_faces; ++i)
      PGs._type[i] = (face_count[i] == 1) ? FLAG : INNER;
  
    return 0;
  }
  
  ///
  E_Int flag_external_phs(E_Int FLAG)
  {
    PHs.updateFacets();
    
    E_Int nb_facets, nb_elts(PHs.size());
    
    // External elements are those with at least one external facets
    PHs._type.resize(nb_elts, INNER);
    for (E_Int i = 0; i < nb_elts; ++i)
    {
      nb_facets = PHs.stride(i);
      for (E_Int j = 0; j < nb_facets; ++j)
      {
        const E_Int& pgi = PHs.get_facet(i,j);
        if (PGs._type[pgi - 1] == FLAG)
        {
          PHs._type[i] = FLAG;
          break;
        }
      }
    }
  
    return 0;
  }
  
  ///
  E_Int remove_pgs(const Vector_t<E_Int>& toremove)
  {
    Vector_t<E_Int> pgnids, phnids;
    return remove_pgs(toremove, pgnids, phnids);
  }
  ///
  E_Int remove_pgs(const Vector_t<E_Int>& toremove, Vector_t<E_Int>& pgnids, Vector_t<E_Int>& phnids)
  {
    // WARNING : index start dependant : 0
    if (toremove.empty())
      return 0;
    //
    /*std::cout << "stats after remove_entities" << std::endl;
    std::cout << "nb pgs : " << PGs.size() << std::endl;
    std::cout <<"pg max " << PGs.get_facets_max_id() << std::endl;
    std::cout << "nb phs : " << PHs.size() << std::endl;
    std::cout <<"ph range " << PHs.get_facets_max_id() << std::endl;
    std::cout << "//////////////////////////////////" << std::endl;*/
        
    pgnids.clear();
    PGs.remove_entities(toremove, pgnids); //0-based nids
    PGs.updateFacets();
    
    /*std::cout << "stats after remove_entities" << std::endl;
    std::cout << "nb pgs : " << PGs.size() << std::endl;
    std::cout <<"pg max " << PGs.get_facets_max_id() << std::endl;
    std::cout << "nb phs : " << PHs.size() << std::endl;
    std::cout <<"ph range " << PHs.get_facets_max_id() << std::endl;
    std::cout << "//////////////////////////////////" << std::endl;*/
    
    //std::cout << "nids size " << nids.size() << std::endl;

    PHs.remove_facets(pgnids, phnids); //0-based nids
    
    //std::cout << "remove_pgs : 3 " << std::endl;
    return 0;
  }
  
  ///
  E_Int build_ph_neighborhood(ngon_unit& neighbor, std::vector<bool>& wall)
  {
    neighbor.clear();
    PHs.updateFacets();
    neighbor=PHs; //same molecules as the connectivity
    neighbor.reset_facets();
    E_Int maxPG = PHs.get_facets_max_id();
    
    Vector_t<E_Int> neigh(maxPG, E_IDX_NONE);
    E_Int nb_elts(PHs.size()), nb_facets, fid;
    E_Int count(0); // 2-pass required to fill both sides of the neighbor info
    while (count++ != 2)
    {
      for (E_Int eid = 0; eid < nb_elts; ++eid)
      {
        nb_facets = PHs.stride(eid);
        for (E_Int j = 0; j < nb_facets; ++j)
        {
          fid = PHs.get_facet(eid,j)-1;
          if (neigh[fid] != E_IDX_NONE) //2nd pass
            neighbor.get_facet(eid,j) = neigh[fid];
          if (!wall[fid]) neigh[fid] = eid;
        }
      }
    }
    return 0;
  }
  
  ///
  E_Int build_ph_neighborhood(ngon_unit& neighbor) const
  {
    neighbor.clear();
    
    if (PHs._NGON.empty())
      return 1;
    
    PHs.updateFacets();
    neighbor=PHs; //same molecules as the connectivity
    neighbor.reset_facets();
    E_Int maxPG = PHs.get_facets_max_id();
    if (maxPG < 0) return 1;

    Vector_t<E_Int> neigh(maxPG, E_IDX_NONE);
    E_Int nb_elts(PHs.size()), nb_facets, fid;
    E_Int count(0); // 2-pass required to fill both sides of the neighbor info
    while (count++ != 2)
    {
      for (E_Int eid = 0; eid < nb_elts; ++eid)
      {
        nb_facets = PHs.stride(eid);
        for (E_Int j = 0; j < nb_facets; ++j)
        {
          fid = PHs.get_facet(eid, j)-1;
          E_Int& eidn = neigh[fid];
          E_Int& Kn = neighbor.get_facet(eid, j);
          if (eidn != E_IDX_NONE && eidn != eid)
            Kn = eidn;
          neigh[fid] = eid;
        }
      }
    }
    return 0;
  }

  ///
  static E_Int build_pg_neighborhood(const ngon_unit& PGs, ngon_unit& neighbor, const std::set<K_MESH::NO_Edge>* wall=0/*extra specified walls*/)
  {
    neighbor.clear();
    PGs.updateFacets();
    neighbor = PGs; //same molecules as the connectivity
    neighbor.reset_facets();

    std::map<K_MESH::NO_Edge, std::pair<K_MESH::Edge, K_MESH::Edge> > noe_to_bound;//oriented edge to its position in the PG
    K_MESH::NO_Edge E;
    E_Int nb_pgs(PGs.size()), nb_facets;
    const E_Int* pN;
    std::map<K_MESH::NO_Edge, std::pair<K_MESH::Edge, K_MESH::Edge> >::iterator it;
    for (E_Int pid = 0; pid < nb_pgs; ++pid)
    {
      nb_facets = PGs.stride(pid);
      pN = PGs.get_facets_ptr(pid);
      for (E_Int j = 0; j < nb_facets; ++j)
      {
        E.setNodes(*(pN + j), *(pN + (j + 1) % nb_facets));
        it = noe_to_bound.find(E);
        if (it == noe_to_bound.end())
        {
          //it = noe_to_bound[E];
          noe_to_bound[E].first = K_MESH::Edge(pid, j);
          noe_to_bound[E].second = K_MESH::Edge(E_IDX_NONE, E_IDX_NONE);
        }
        else
          it->second.second = K_MESH::Edge(pid, j);
      }
    }

    if (wall)
      for (std::set<K_MESH::NO_Edge>::const_iterator i = wall->begin(); i != wall->end(); ++i) noe_to_bound.erase(*i);

    std::map<K_MESH::NO_Edge, std::pair<K_MESH::Edge, K_MESH::Edge> >::iterator itE = noe_to_bound.end();
    for (it= noe_to_bound.begin(); it != itE; ++it)
    {
      if (it->second.second.node(0) == E_IDX_NONE)
        continue;

      const E_Int& pg1 = it->second.first.node(0);
      const E_Int& n1 = it->second.first.node(1);
      const E_Int& pg2 = it->second.second.node(0);
      const E_Int& n2 = it->second.second.node(1);

      neighbor.get_facet(pg1, n1) = pg2;
      neighbor.get_facet(pg2, n2) = pg1;
    }

    return 0;
  }
  
  ///
  static void refine_pgs(const std::map<K_MESH::NO_Edge, Vector_t<E_Int> >& edge_to_refined_edge, ngon_unit& PGs)
  {
    ngon_unit refinedPGs;
    Vector_t<E_Int> pg_molec;
    
    for (E_Int i = 0; i < PGs.size(); ++i)
    {
      refine_pg(PGs.get_facets_ptr(i), PGs.stride(i), edge_to_refined_edge, pg_molec);
      refinedPGs.add(pg_molec);
    }
    refinedPGs._type = PGs._type;  //work around to not loose them
    refinedPGs._ancEs = PGs._ancEs;//work around to not loose them
    PGs = refinedPGs; //dirty
  }
  
  ///
  static void refine_pg
  (const E_Int* pg, E_Int sz, const std::map<K_MESH::NO_Edge, Vector_t<E_Int> >& edge_to_refined_edge, Vector_t<E_Int>& pg_molec)
  {
    E_Int Ni, Nj;
    K_MESH::NO_Edge Ei;
    std::map<K_MESH::NO_Edge, Vector_t<E_Int> >::const_iterator it;
    
    pg_molec.clear();
    pg_molec.resize(1, 0);
    for (E_Int i = 0; i < sz; ++i)
    {
      Ni = *(pg+i);
      Nj = *(pg + (i+1)%sz);
      Ei.setNodes(Ni, Nj);
      
      it = edge_to_refined_edge.find(Ei);
      if (it != edge_to_refined_edge.end())
      {
        const Vector_t<E_Int>& nodes = it->second;
        if (Ni == Ei.node(0)) //same orientation
          for (size_t j=0; j< nodes.size()-1; ++j)
            pg_molec.push_back(nodes[j]);
        else
          for (E_Int j=nodes.size()-1; j>0; --j)
            pg_molec.push_back(nodes[j]);
      }
      else
      {
        pg_molec.push_back(Ni);
        //pg_molec.push_back(Nj);
      }
    }
    
    pg_molec[0]=pg_molec.size()-1;
  }
  
  static bool is_closed(const ngon_t& ngio, E_Int PHi, std::set<K_MESH::NO_Edge>& free_edges)
  {
    //std::cout << ngio.PHs << std::endl;
    E_Int nb_pgs = ngio.PHs.stride(PHi);
    const E_Int* const ptrPG = ngio.PHs.get_facets_ptr(PHi);
    E_Int PGi, nb_nodes;
    K_MESH::NO_Edge noE;
    const  E_Int *pNi;
    
    free_edges.clear();
    
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      PGi = *(ptrPG+i)-1;
      nb_nodes=ngio.PGs.stride(PGi);
      pNi = ngio.PGs.get_facets_ptr(PGi);
      //
      for (E_Int i=0; i < nb_nodes; ++i)
      {
        const E_Int& Ni = *(pNi+i);
        const E_Int& Nj = *(pNi+(i+1)%nb_nodes);
             
        noE.setNodes(Ni, Nj);
        
        if (free_edges.find(noE) != free_edges.end()) //if already in, closed border
          free_edges.erase(noE);
        else
          free_edges.insert(noE);
      }
    }
    return free_edges.empty();
  }
  
  /// Close a PH having free edges (GEOM algo).
  static void close_phs(ngon_t& ngio, const K_FLD::FloatArray& coord)
  {
    //WARNING : algo valid only for quasi-straight BAR openings (because the poinst are sorted by distance).
    // TYPICAL USE : after Cassiopee NGON conformization of an octree (new points on a subdivided face are not propagated to the attached faces..)
    
    std::map<K_MESH::NO_Edge, Vector_t<E_Int> > edge_to_refined_edge;
    ngon_unit refinedPGs;
    Vector_t<E_Int> pg_molec, refined_edge;
    std::set<K_MESH::NO_Edge> edges;
    K_MESH::NO_Edge noE;
    K_FLD::IntArray conn;
    
    ngio.PGs.updateFacets();
    ngio.PHs.updateFacets();
    
    // Build the refined edges       
    E_Int nb_phs(ngio.PHs.size());
    //E_Int PGi;
    for (E_Int PHi = 0; PHi < nb_phs; ++PHi)
    {      
      //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PH("/home/slandier/tmp/open_borders.tp", coord, ngio, PHi);
      //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::highlight_PH(ngio, coord, PHi);
      
      // Get the edges to refine (for the entire PH).
      //E_Int nb_pgs = ngio.PHs.stride(PHi);
      //const E_Int* const ptrPG = ngio.PHs.get_facets_ptr(PHi);
      
      if (is_closed(ngio, PHi, edges)) continue;
      
      /*if (edges.empty()) // PH is closed, nothing to do.
      {
        for (E_Int i=0; i < nb_pgs; ++i)
        {
          PGi = *(ptrPG+i)-1;
          refinedPGs.__add(ngio.PGs, PGi);
        }
        continue;
      }*/
      
#ifdef DEBUG_NGON_T
      /*{
        K_FLD::IntArray connectE;
        std::set<K_MESH::NO_Edge>::const_iterator it;
        for (it=edges.begin(); it!=edges.end(); ++it)
        {
          K_MESH::Edge e(it->node(0), it->node(1));
          connectE.pushBack(e.begin(), e.end());
        }
        connectE.shift(-1);  
        MIO::write("/home/slandier/tmp/free_edges.mesh", coord, connectE, "BAR");
      }*/
#endif
      
      //
      conn.clear();
      conn.resize(2, edges.size());
      E_Int i(0);
      for (std::set<K_MESH::NO_Edge>::const_iterator ii=edges.begin(); ii!= edges.end(); ++ii, ++i)
      { conn(0,i)=ii->node(0)-1; conn(1,i)=ii->node(1)-1; }
      
      // Split into single contours
      typedef std::map<E_Int, std::vector<E_Int> > ntn_t;
      ntn_t node_to_nodes;
      BARSplitter::get_node_to_nodes(conn, node_to_nodes);
      
      // For each contour
      for (ntn_t::const_iterator n=node_to_nodes.begin(); n!= node_to_nodes.end(); ++n)
      {
        if (n->second.size() != 2) // not a splitting node
          continue;

        const E_Int& Nsplit = n->first ;

        refined_edge.clear();
        noE.setNodes(1+n->second[0], 1+n->second[1]);
        refined_edge.push_back(noE.node(0));
        refined_edge.push_back(Nsplit+1);
        refined_edge.push_back(noE.node(1));
        edge_to_refined_edge[noE]=refined_edge;
      } 
    }
    
    // Refine the PGs ( <=> close the PHs )
    E_Int nb_pgs = ngio.PGs.size();
    for (E_Int PGi = 0; PGi < nb_pgs; ++PGi)
    {
      refine_pg(ngio.PGs.get_facets_ptr(PGi), ngio.PGs.stride(PGi), edge_to_refined_edge, pg_molec);
      refinedPGs.add(pg_molec);
    }
    
    refinedPGs._type = ngio.PGs._type;  // hack to preserve flags (externality)
    refinedPGs._ancEs = ngio.PGs._ancEs;// hack
    ngio.PGs = refinedPGs;
    ngio.PGs.updateFacets();//fixme : required ?
    
  }
  
  /// Change the node indice to reference the same when duplicated exist
  E_Int join_phs(const K_FLD::FloatArray& coord, E_Float tolerance = E_EPSILON)
  {
    if (PGs.size() == 0) return 0;
    
    PGs.updateFacets();
    
    K_FLD::ArrayAccessor<K_FLD::FloatArray> ca(coord);
    Vector_t<E_Int> nids;
    E_Int nb_merges = ::merge(ca, tolerance, nids);
    if (nb_merges)
      PGs.change_indices(nids, true/*zero based*/);
    return nb_merges;
  }
  
  /// Change the node indice to reference the same when duplicated exist (FldArrayF)
  void join_phs(const K_FLD::FldArrayF& coord, E_Int px, E_Int py, E_Int pz, E_Float tolerance = E_EPSILON)
  {
    if (PGs.size() == 0)
      return;
    
    PGs.updateFacets();
    
    K_FLD::ArrayAccessor<K_FLD::FldArrayF> ca(coord, px, py, pz);
    Vector_t<E_Int> nids;
    E_Int nb_merges = ::merge(ca, tolerance, nids);
    if (nb_merges)
      PGs.change_indices(nids, true/*zero based*/);
  }
  /*   
  ///
  E_Int unique_nodes_ph(E_Int PHi, Vector_t<E_Int>& indices, bool zerobased=false) const
  {
    PHs.updateFacets();
    PGs.updateFacets();
    
    indices.clear();
    
    std::set<E_Int> tmp;
    
    E_Int nb_pgs(PHs.stride(PHi)), nb_nodes;
   
    for (E_Int j = 0; j < nb_pgs; ++j)
    {
      const E_Int& pgi = PHs.get_facet(PHi,j);
      nb_nodes = PGs.stride(pgi-1);
      
      for (E_Int k = 0; k < nb_nodes; ++k)
        tmp.insert(PGs.get_facet(pgi-1,k));
    }
    
    indices.insert(indices.end(), tmp.begin(), tmp.end());
    if (zerobased)
      ngon_unit::shift(indices, -1);
    
    return 0;
  }*/

  E_Int nodes_ph(E_Int PHi, Vector_t<E_Int>& indices, bool zerobased = false) const
  {
    PHs.updateFacets();
    PGs.updateFacets();

    indices.clear();

    E_Int nb_pgs(PHs.stride(PHi)), nb_nodes;

    for (E_Int j = 0; j < nb_pgs; ++j)
    {
      const E_Int& pgi = PHs.get_facet(PHi, j);
      nb_nodes = PGs.stride(pgi - 1);

      for (E_Int k = 0; k < nb_nodes; ++k)
        indices.push_back(PGs.get_facet(pgi - 1, k));
    }

    if (zerobased)
      ngon_unit::shift(indices, -1);

    return 0;
  }

  ///
  //E_Int unique_nodes_phs(Vector_t<E_Int>& indices, bool zerobased=false) const
  //{
  //  PHs.updateFacets();
  //  indices.clear();
  //  
  //  Vector_t<E_Int> inds;
  //  std::set<E_Int> tmp;
  //  
  //  E_Int nb_facets, nb_elts(PHs.size());
  //  
  //  for (E_Int i = 0; i < nb_elts; ++i)
  //  {
  //    nodes_ph(i, inds, false); //false to do the shift only once
  //    tmp.insert(inds.begin(), inds.end());
  //  }
  //  
  //  indices.insert(indices.end(), tmp.begin(), tmp.end());
  //  if (zerobased)
  //    ngon_unit::shift(indices, -1);
  //  
  //  return 0;
  //}
  
  ///
  E_Int unique_nodes_ph(E_Int PHi, Vector_t<E_Int>& indices, bool zerobased=false) const
  {
    PHs.updateFacets();
    indices.clear();
    
    Vector_t<E_Int> inds;
    std::set<E_Int> tmp;
    
    E_Int nb_facets;
    nodes_ph(PHi, inds, false); //false to do the shift only once
    tmp.insert(inds.begin(), inds.end());
   
    indices.insert(indices.end(), tmp.begin(), tmp.end());
    if (zerobased)
      ngon_unit::shift(indices, -1);
    
    return 0;
  }
  
  /// triangulate specified PGs
  template <typename TriangulatorType>
  static E_Int triangulate_pgs
    (const ngon_unit& PGs, const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, Vector_t<E_Int>& colors, const Vector_t<bool>* process=0) 
  {
    connectT3.clear();
    colors.clear();
    
    PGs.updateFacets();
    
#ifdef DEBUG_NGON_T
    ngon_unit pg_err, pg_ok;
#endif
    
    E_Int sz = PGs.size(), err(0);
    TriangulatorType dt;

    for (E_Int i=0; (i<sz) && !err ; ++i)
    {
      if (process && (*process)[i] == false)
        continue;

      err = K_MESH::Polygon::triangulate(dt, coord, PGs.get_facets_ptr(i), PGs.stride(i), 1/*index start*/, connectT3);
      
      colors.resize(connectT3.getSize(), i);
      
#ifdef DEBUG_NGON_T
      if (err)
        pg_err.__add(PGs, i);
      else
        pg_ok.__add(PGs, i);
#endif
    }
  
#ifdef DEBUG_NGON_T
  if (pg_err.size() !=0)
  {
    const char* fname = "err.mesh";
    pg_err.updateFacets();
    ngon_t tmp(pg_err, true);
    K_FLD::IntArray conn;
    tmp.export_to_array(conn);
    MIO::write("pg_err.tp", coord, conn, "NGON");
  }
  if (pg_ok.size() !=0)
  {
    const char* fname = "ok.mesh";
    pg_ok.updateFacets();
    ngon_t tmp(pg_ok, true);
    K_FLD::IntArray conn;
    tmp.export_to_array(conn);
    MIO::write("pg_ok.tp", coord, conn, "NGON");
  }
#endif
    return err;
  }

  /// triangulate specified PGs
  template <typename TriangulatorType>
  static E_Int triangulate_pgs
  (const ngon_unit& PGs, Vector_t<E_Int>& PGis, E_Int index_start, const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3)
  {
    connectT3.clear();
    Vector_t<E_Int> colors; //fixme : hpc
    TriangulatorType dt;
    PGs.updateFacets();
    
#ifdef DEBUG_NGON_T
    ngon_unit pg_err, pg_ok;
#endif
    
    E_Int sz = PGis.size(), err(0);
    
    for (E_Int i=0; (i<sz) && !err ; ++i)
    {
      E_Int PGi = PGis[i] - index_start;
      err = K_MESH::Polygon::triangulate<TriangulatorType>(dt, coord, PGs.get_facets_ptr(PGi), PGs.stride(PGi), 1/*index start*/, connectT3);
      
#ifdef DEBUG_NGON_T
      if (err)
        pg_err.__add(PGs, PGi);
      else
        pg_ok.__add(PGs, PGi);
#endif
    }
  
#ifdef DEBUG_NGON_T
  if (pg_err.size() !=0)
  {
    const char* fname = "err.mesh";
    pg_err.updateFacets();
    ngon_t tmp(pg_err, true);
    K_FLD::IntArray conn;
    tmp.export_to_array(conn);
    MIO::write("pg_err.tp", coord, conn, "NGON");
  }
  if (pg_ok.size() !=0)
  {
    const char* fname = "ok.mesh";
    pg_ok.updateFacets();
    ngon_t tmp(pg_ok, true);
    K_FLD::IntArray conn;
    tmp.export_to_array(conn);
    MIO::write("pg_ok.tp", coord, conn, "NGON");
  }
#endif
    return err;
  }

  template<typename TriangulatorType>
  static E_Int triangulate_pgs
    (const ngon_unit& PGs, const K_FLD::FloatArray& crd, ngon_unit& tri_pgs, Vector_t<E_Int>& oids, E_Float dummy_tol = 0., const Vector_t<bool>* to_process = 0)
  {
    tri_pgs.clear();
    oids.clear();

    // triangulate specified PGs
    TriangulatorType dt;
    K_FLD::IntArray connectT3;
    E_Int err = ngon_t::triangulate_pgs<TriangulatorType>(PGs, crd, connectT3, oids, to_process);
    // convert triangulation to ngon_unit
    ngon_unit::convert_fixed_stride_to_ngon_unit(connectT3, tri_pgs);
    tri_pgs.shift(1);
    //update history
    if (PGs._ancEs.cols())
    {
      for (E_Int i = 0; i < oids.size(); ++i)
        tri_pgs._ancEs.pushBack(PGs._ancEs.col(oids[i]), PGs._ancEs.col(oids[i]) + 2);
    }
    if (!PGs._type.empty())
    {
      for (E_Int i = 0; i < oids.size(); ++i)
        tri_pgs._type.push_back(PGs._type[oids[i]]);
    }
    return err;
  }

  

  /// transform specified PGs and synchronize the PHs
  template <typename TriangulatorType>
  static E_Int transform_pgs
    (const K_FLD::FloatArray& crd, const ngon_t& ngi, Vector_t<bool>&to_process, ngon_t& ngo, pg_transfo_func transFunc, E_Float optional_tol/*for convexify only*/)
  {
    ngon_unit split_pgs/*list of split pgs*/, pgs_split_graph/*split ids per initial pg*/;

    Vector_t<E_Int> oids;
    transFunc(ngi.PGs, crd, split_pgs, oids, optional_tol, &to_process);
    
    // create the split graph
    E_Int nb_ini_pgs = ngi.PGs.size();
    ngi.PGs.create_pgs_changes_graph(oids, true, nb_ini_pgs/*shift*/,pgs_split_graph);// shift in order to be synch with appending split_pgs
    
    // Apply it
    ngo = ngi;
    ngo.apply_pgs_changes_to_PHs(pgs_split_graph);
    ngo.PGs.append(split_pgs);//hack : must be done after apply_pgs_changes_to_PHs
    Vector_t<E_Int> pgnids, phnids;
    ngo.remove_unreferenced_pgs(pgnids, phnids);
    return 0;
  }

  void apply_pgs_changes_to_PHs(const ngon_unit& pgs_split_graph)
  {
    // pgs_split_graph is sized as PGs, if nothing to change for a given PG, molec is (1,E_IDX_NONE), if to remove (0) 
    assert(PGs.size() == pgs_split_graph.size());

    ngon_unit new_phs;
    Vector_t<E_Int> molecPH;

    pgs_split_graph.updateFacets();
    PHs.updateFacets();
    bool transfer_ancs = (PHs._ancEs.cols() != 0);
    bool transfer_type = (!PHs._type.empty());

    E_Int nb_phs = PHs.size();

    if (transfer_ancs) new_phs._ancEs.reserve(2, nb_phs);
    if (transfer_type) new_phs._type.reserve(nb_phs);

    for (E_Int i = 0; i < nb_phs; ++i)
    {
      const E_Int* pPGi = PHs.get_facets_ptr(i);
      E_Int nb_pgs = PHs.stride(i);

      molecPH.clear();

      for (E_Int p = 0; p < nb_pgs; ++p)
      {
        E_Int PGi = *(pPGi + p) - 1;
        E_Int nb_split = pgs_split_graph.stride(PGi);
        const E_Int* new_facets = pgs_split_graph.get_facets_ptr(PGi);

        if (nb_split == 0) //to remove
          continue;
        if (nb_split == 1 && *new_facets == E_IDX_NONE) //keep it as unchanged
          molecPH.push_back(PGi + 1);
        else
          molecPH.insert(molecPH.end(), new_facets, new_facets + nb_split);
      }

      E_Int sz = molecPH.size();
      if (sz == 0) // PH is gone
        continue;

      new_phs.add(molecPH.size(), &molecPH[0]);

      if (transfer_ancs)
        new_phs._ancEs.pushBack(PHs._ancEs.col(i), PHs._ancEs.col(i) + 2);
      if (transfer_type)
        new_phs._type.push_back(PHs._type[i]);
    }

    PHs = new_phs;
    PHs.updateFacets();
  }

  ///
  template <typename TriangulatorType>
  static E_Int triangulate_external_pgs(K_FLD::FloatArray& crd, ngon_t& ngi, ngon_t& ngo)
  {
    //detect exterior faces
    ngi.flag_external_pgs(INITIAL_SKIN);

    E_Int nb_pgs = ngi.PGs.size();
    Vector_t<bool> toprocess(nb_pgs, false);
    for (size_t i = 0; i < nb_pgs; ++i)
      toprocess[i] = (ngi.PGs._type[i] == INITIAL_SKIN);

    // triangulate any specified PG
    ngon_t::transform_pgs<TriangulatorType>(crd, ngi, toprocess, ngo, triangulate_pgs<TriangulatorType>, 0./*dummy*/);
  
    return 0;
  }

  ///
  template <typename TriangulatorType>
  static E_Int triangulate_ph
  (const TriangulatorType& triangulator, const ngon_t& ng, E_Int PHi, const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3)
  {
    connectT3.clear();

    ng.PHs.updateFacets();
    ng.PGs.updateFacets();
    
    const E_Int* begin = ng.PHs.get_facets_ptr(PHi);
    E_Int nb_faces = ng.PHs.stride(PHi);
    
    Vector_t<E_Int> PGis;
    PGis.insert(PGis.end(), begin, begin + nb_faces);

    return triangulate_pgs<TriangulatorType>(ng.PGs, PGis, 1 /*index_start*/, coord, connectT3);
  }

  ///
  template <typename Coordinate_t>
  static void compact_to_used_nodes (ngon_unit& PGs, Coordinate_t& coord)
  {
    Vector_t<E_Int> unic_nodes;
    PGs.unique_indices(unic_nodes);
  
    Vector_t<bool> flag(coord.getSize(), false);
  
    size_t sz(unic_nodes.size());
    for (size_t i = 0; i < sz; ++i) flag[unic_nodes[i]-1]=true; // indices start at 1.
  
    Vector_t<E_Int> newIDs;
    E_Int nb_removed = Coordinate_t::compact(coord, flag, newIDs);
    if (nb_removed == 0) return;
  
    PGs.change_indices(newIDs, true/*zero based*/); 
  }
  
  ///fixme
  template <typename Coordinate_t>
  static void compress_to_used_nodes (ngon_unit& PGs, Coordinate_t& coord)
  {
    Vector_t<E_Int> unic_nodes;
    PGs.unique_indices(unic_nodes);
     
    Vector_t<E_Int> nids(coord.getSize(), E_IDX_NONE), nidscpy;
  
    size_t sz(unic_nodes.size());
    for (size_t i = 0; i < sz; ++i) nids[unic_nodes[i]-1]=unic_nodes[i]-1; // indices start at 1.
    
    nidscpy=nids;
    K_CONNECT::valid pred(nidscpy);
    
    E_Int nb_removed = K_CONNECT::IdTool::compress(nids, pred);
    if (nb_removed == 0)
      return;
    
    nb_removed = K_CONNECT::IdTool::compress(coord, pred);
  
    PGs.change_indices(nids, true/*zero based*/); 
  }
   
  ///
  template <typename CoordAccType>
  bool remove_duplicated_pgs (const CoordAccType& coord)
  {
    bool found = replace_duplicated_pgs(coord); //detect the matching and replace the ids with the first one found
    if (found)
      PHs.remove_duplicated(); //clean
    return found;
  }
  
  bool remove_duplicated_edges()
  {
    bool found = replace_duplicated_edges();
    if (found)
      PHs.remove_duplicated();//clean
    return found;
  }
  
  bool remove_duplicated_nodes()
  {
    bool found = replace_duplicated_nodes();
    if (found)
      PHs.remove_duplicated();//clean
    return found;
  }
  
  ///
  template <typename CoordAccType>
  bool replace_duplicated_pgs (const CoordAccType& coord)
  {
    if (PHs.size()*PGs.size() == 0)
      return false;
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int nb_pgs, nb_phs(PHs.size()), pgi, count(0);
    E_Float G[3];
    K_FLD::FloatArray barys;
    barys.reserve(3, PGs.size());
        
    Vector_t<E_Int> ISO(PHs.get_facets_max_id(), E_IDX_NONE), PGI;
    
    // Loop over PH rather than directly over PGs to consider only used PGs
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      nb_pgs = PHs.stride(i);
      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        pgi = PHs.get_facet(i, j)-1; 
        if (ISO[pgi] != E_IDX_NONE) //already stored
          continue;
        E_Int s = PGs.stride(pgi);
        K_MESH::Polygon PGi(PGs.get_facets_ptr(pgi), s, -1/*to make it zero based*/);
        PGi.iso_barycenter<CoordAccType, 3>(coord, G);
        barys.pushBack(G, G+3);
        
        ISO[pgi]=count++;
      }
    }
    
    K_CONNECT::IdTool::reverse_indirection(ISO, PGI);
    
    // Merge
    Vector_t<E_Int> nids;
    K_FLD::ArrayAccessor<K_FLD::FloatArray> cab(barys);
    E_Int nmerges = ::merge(cab, E_EPSILON, nids);
    
    if (!nmerges) //no duplicated faces.
      return false;
    
    // check if matching isoG means really identical PGs
    bool found = false;
    for (size_t i=0; i < nids.size(); ++i)
    {
      if (nids[i] == (E_Int)i)
        continue;
      
      const E_Int &pgi1 = PGI[i];
      const E_Int &pgi2 = PGI[nids[i]];
      
      if (PGs.stride(pgi1) != PGs.stride(pgi2)) // fast treatment : not the same number of nodes
      {
        nids[i]=i;
        continue;
      }
      if (!K_CONNECT::IdTool::equal(PGs.get_facets_ptr(pgi1), PGs.get_facets_ptr(pgi2), PGs.stride(pgi1), true/*permut*/, false/*strict orient*/))
        nids[i]=i;
      else
        found=true;
    }
    
    if (!found)
      return false;
    
    K_CONNECT::IdTool::propagate(nids, ISO);
    
    
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      nb_pgs = PHs.stride(i);
      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        E_Int & pgi = PHs.get_facet(i, j); 
        pgi = PGI[ISO[pgi-1]]+1;
      }
    }
    return true;
  }
  
  /// surfacic version of replace_duplicated_pgs : here PHs are polygones and PGs are edges..
  bool replace_duplicated_edges()
  {
    if (PHs.size()*PGs.size() == 0)
      return false;
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int nb_edges(PGs.size()); // nb_polys(PHs.size());
    std::map<K_MESH::NO_Edge, E_Int> edge_to_id;
    Vector_t<E_Int> nids(nb_edges, E_IDX_NONE);
    K_MESH::NO_Edge E;
    bool found = false;
    for (E_Int i = 0; i < nb_edges; ++i)
    {
      const E_Int* ptr = PGs.get_facets_ptr(i);
      E.setNodes(ptr);
      if (edge_to_id.find(E) != edge_to_id.end())
      {
        nids[i] = edge_to_id[E];
        found = true;
      }
      else
      {
        edge_to_id[E]=i;
        nids[i]=i;
      }
    }
    
    if (!found)
      return false;
    
    PHs.change_indices(nids, true/*zerobased*/);
    return true;
  }
  
  /// lineic version of replace_duplicated_pgs : here PHs are edges and PGs are nodes..
  bool replace_duplicated_nodes()
  {
    if (PHs.size()*PGs.size() == 0)
      return false;
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int Ni, nb_nodes(PGs.size());
    Vector_t<E_Int> nids(nb_nodes, E_IDX_NONE);
    Vector_t<E_Int> pgid(PGs.get_facets_max_id(), E_IDX_NONE);
    bool found = false;
    for (E_Int i = 0; i < nb_nodes; ++i)
    {
      const E_Int* ptr = PGs.get_facets_ptr(i);
      Ni = *ptr-1;
      if (pgid[Ni] != E_IDX_NONE)
      {
        nids[i] = pgid[Ni];
        found = true;
      }
      else
        nids[i]=pgid[Ni]=i;
    }
    
    if (!found)
      return false;
    
    PHs.change_indices(nids, true/*zerobased*/);
    return true;
  }
  
  ///
  template <typename CoordAccType>
  E_Int detect_phs_with_same_centroid (const CoordAccType& coord, Vector_t<E_Int>& nids)
  {
    if (PHs.size()*PGs.size() == 0)
      return false;
    
    PHs.updateFacets();
    PGs.updateFacets();
    
    E_Int nb_pgs, nb_phs(PHs.size()), pgi, count(0);
    E_Float g[3], G[3], K;
    K_FLD::FloatArray barys;
    barys.reserve(3, nb_phs);
    
    //
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      nb_pgs = PHs.stride(i);
      G[0]=G[1]=G[2]=0.;
      for (E_Int j = 0; j < nb_pgs; ++j)
      {
        pgi = PHs.get_facet(i, j)-1; 
        E_Int s = PGs.stride(pgi);
        K_MESH::Polygon PGi(PGs.get_facets_ptr(pgi), s, -1/*to make it zero based*/);
        PGi.iso_barycenter<CoordAccType, 3>(coord, g);
        
        for (size_t k=0; k < 3; ++k)
          G[k] += g[k];
      }
      K = 1./(E_Float)nb_pgs;
      for (size_t k=0; k < 3; ++k)
        G[k] *= K;
        
      barys.pushBack(G, G+3);
    }
    
    // Merge
    nids.clear();
    K_FLD::ArrayAccessor<K_FLD::FloatArray> cab(barys);
 
    return ::merge(cab, E_EPSILON, nids);
  }
  
  ///
  E_Int remove_degenerated_pgs(Vector_t<E_Int>& pgnids, Vector_t<E_Int>& phnids)
  {
    Vector_t<E_Int> pgindices;
    PGs.get_degenerated(pgindices); // Degenrated edge(<2) or Polygons (<3)
    remove_pgs(pgindices, pgnids, phnids); //sync the PHs accordingly
    PGs.updateFacets();//fixme : required ?
    return pgindices.size();
  }
  
  ///
  E_Int remove_unreferenced_pgs(Vector_t<E_Int>& pgnids, Vector_t<E_Int>& phnids)
  {
    PHs.updateFacets();
    PGs.updateFacets();

#ifdef DEBUG_NGON_T
    assert(is_consistent());
#endif
    
    Vector_t<E_Int> visited(PGs.size(),0);
    for (E_Int i=0; i < PHs.size(); ++i)
      for (E_Int j=0; j < PHs.stride(i); ++j)
        ++visited[PHs.get_facet(i,j)-1];
   
    Vector_t<E_Int> unused;
    for (size_t i=0; i < visited.size(); ++i)if (!visited[i])unused.push_back(i);
    remove_pgs(unused, pgnids, phnids);
    return (E_Int)unused.size();
  }
  
  // "PRIVATE" METHODS
  ///
  static void __get_selected_ids(const ngon_unit& ngon_in, const Vector_t<bool>& keep,
                                 Vector_t<E_Int>& new_ids, Vector_t<E_Int>& stack_ids)
 {
    size_t sz0 = ngon_in.size();
    assert(sz0 == keep.size());
  
    ngon_in.updateFacets();
  
    size_t old_facets_sz = ngon_in.get_facets_max_id();
    
    new_ids.resize(old_facets_sz+1, E_IDX_NONE);
    stack_ids.clear();
  
    size_t count(0);
    E_Int pos, sz, Fij;
  
    for (size_t i = 0; i < sz0; ++i)
    {
      if (!keep[i])
        continue;
    
      pos = ngon_in._facet[i];
      sz = ngon_in._NGON[pos];
    
      for (E_Int j = 0; j < sz; ++j)
      {
        Fij = ngon_in.get_facet(i,j);
        if (new_ids[Fij] == E_IDX_NONE)
        {
          new_ids[Fij]=++count;
          stack_ids.push_back(Fij);
        }
      }
    }
  }

  ///
  static void __create_selected_phs
  (const ngon_unit& ngon_in, const Vector_t<bool>& keep,
   const Vector_t<E_Int>& nids, ngon_unit& ngon_out)
  {
    size_t ksz(keep.size());
    E_Int pos, sz, Fij;
  
    ngon_out.clear();
    ngon_out._NGON.resize(2, 0);
  
    ngon_in.updateFacets();//refresh in case
  
    for (size_t i = 0; i < ksz; ++i)
    {
      if (!keep[i]) continue;
    
      ++ngon_out._NGON[0];
    
      pos = ngon_in._facet[i];
      sz = ngon_in._NGON[pos];
      ngon_out._NGON.push_back(sz);

      for (E_Int j = 0; j < sz; ++j)
      {
        Fij = ngon_in.get_facet(i,j);
        ngon_out._NGON.push_back(nids[Fij]);
      } 
    }
  
    ngon_out._NGON[1]=ngon_out._NGON.size()-2;
    
    //transfer externality
    ngon_out.__transfer_attributes(ngon_in, keep);
    
  }

  ///
  static void __create_selected_pgs
  (const ngon_unit& ngon_in, const Vector_t<E_Int>& stack_ids, ngon_unit& ngon_out)
  {
    size_t ksz(stack_ids.size());
    E_Int pos,sz, Fij;
  
    ksz = stack_ids.size();
  
    ngon_out.clear(); 
    ngon_out._NGON.resize(2, 0);
    ngon_out._NGON[0]=ksz;
  
    ngon_in.updateFacets();//refresh in case
  
    bool tranfer_type = !ngon_in._type.empty();
    bool transfer_anc = ngon_in._ancEs.cols();
  
    for (size_t i = 0; i < ksz; ++i)
    {
      const E_Int & Fi = stack_ids[i];
    
      pos = ngon_in._facet[Fi-1];
      sz = ngon_in._NGON[pos];
      ngon_out._NGON.push_back(sz);
      //transfer externality
      if (tranfer_type)
        ngon_out._type.push_back(ngon_in._type[Fi-1]);
      if (transfer_anc)
        ngon_out._ancEs.pushBack(ngon_in._ancEs.col(Fi - 1), ngon_in._ancEs.col(Fi - 1) + 2);
    
      for (E_Int j = 0; j < sz; ++j)
      {
        Fij = ngon_in.get_facet(Fi-1,j);
        ngon_out._NGON.push_back(Fij);
      }
    }
  
    ngon_out._NGON[1]=ngon_out._NGON.size()-2;
  }
  
  ///
  E_Int keep_PHs_having_PGs(const Vector_t<E_Int>& PGs_tokeep)
  {
    // index start dependant
    Vector_t<bool> keep(PGs.size(), false);
    for (size_t i = 0; i < PGs_tokeep.size(); ++i) keep[PGs_tokeep[i]]=true;

    keep_PHs_having_PGs(keep);

    return 0;
  }
  
  ///
  E_Int keep_PHs_having_PGs(const Vector_t<bool>& keep)
  {
    // index start dependant
    ngon_unit keptPHs;
    PHs.updateFacets();
    bool kp;
    E_Int sz, *pN;
    for (E_Int i = 0; i < PHs.size(); ++i)
    {
      sz = PHs.stride(i);
      pN = PHs.get_facets_ptr(i);
      kp = false;
      for (E_Int j = 0; (j < sz) && !kp; ++j)
        kp = keep[*(pN + j) - 1];

      if (kp)
        keptPHs.__add(PHs, i);
    }
  
    PHs=keptPHs;
    PHs.updateFacets();

    return 0;
  }
  
  ///
  E_Int flag_PHs_having_PGs(const Vector_t<bool>& keepPG, Vector_t<bool>& keepPH) const 
  {
    // index start dependant
    
    PHs.updateFacets();
    
    E_Int sz(PHs.size()), nb_pgs;
    
    keepPH.clear();
    keepPH.resize(sz, false);
    
    for (E_Int i = 0; i < sz; ++i)
    {
      nb_pgs = PHs.stride(i);
      const E_Int* pN = PHs.get_facets_ptr(i);
 
      for (E_Int j = 0; (j < nb_pgs); ++j)
      {
        if (keepPG[*(pN + j) - 1])
        {
          keepPH[i]=true;
          break;
        }
      }
    }
  
    return 0;
  }
  
  /// Warning : The coordinates are not cleaned, only the connectivity.
  static E_Int clean_connectivity
  (ngon_t& NG, const K_FLD::FloatArray& f, E_Int ngon_dim=-1, E_Float tolerance = E_EPSILON)
  {
    E_Int nb_phs0(NG.PHs.size()), nb_pgs0(NG.PGs.size());
    
    if (nb_pgs0 == 0) // fast return
      return 0;
    
    typedef K_FLD::ArrayAccessor<K_FLD::FloatArray> acrd_t;
    acrd_t fcA(f);
  
    NG.PGs.updateFacets();
    NG.PHs.updateFacets();
  
    // type of NGON: volumic , surfacic or lineic?
    if (ngon_dim != 1 && ngon_dim != 2 && ngon_dim != 3) // not a kwnown value
    {
      for (E_Int i=0; i<nb_pgs0; ++i)
        ngon_dim = std::max(ngon_dim, NG.PGs.stride(i));
      
      if (ngon_dim != 1 && ngon_dim !=2)
        ngon_dim=3;//volumic
    }

    // 1- Referencement unique pour les noeuds confondus par kdtree (les doublons seront supprimes a la fin)
    E_Int nb_merges = NG.join_phs(f, tolerance);
  
    // 2- Elimination des faces degenerees
    Vector_t<E_Int> pgnids, phnids; // required to update the history (PG/PH)
    E_Int nb_degen_faces(0), nb_consec_changes(0);
    if (ngon_dim != 1)
    {
      nb_degen_faces = NG.remove_degenerated_pgs(pgnids, phnids);
      nb_consec_changes = NG.PGs.remove_consecutive_duplicated(); //removing duplicated nodes : compact representation
    }
   
    // 3- Faces confondues : identification et suppression des références.
    bool has_dups = false;
    if (ngon_dim ==3) //volumic
      has_dups = NG.remove_duplicated_pgs(fcA);
    else if (ngon_dim == 2) //surfacic
      has_dups = NG.remove_duplicated_edges();
    else // lineic
      has_dups = NG.remove_duplicated_nodes();

    // remove duplicated references to PGs within each elements
    E_Int nb_phs_dups = NG.PHs.remove_duplicated();

    // 4- Elimination des elts degeneres
    Vector_t<E_Int> toremove;
    NG.PHs.get_degenerated(toremove);
    E_Int nb_degen_phs = NG.PHs.remove_entities(toremove, phnids);
    
    // 5- Elimination des elts doubles
    //
  
    // 6- Suppression des faces non referencees
    E_Int nb_unrefs = NG.remove_unreferenced_pgs(pgnids, phnids); //Important, sinon tecplot en ASCII (.tp) n'aime pas.
    
    // 7- Compression du ngon aux seuls noeuds utilises
    //ngon_t::compact_to_used_nodes(NG.PGs, f);


    E_Int nb_phs1 = NG.PHs.size();
    E_Int nb_pgs1 = NG.PGs.size();
    
#ifdef DEBUG_NGON_T
  assert (NG.is_consistent());
#endif

    
    return nb_merges + (nb_phs0-nb_phs1) + (nb_pgs0-nb_pgs1);
  }

  static void edge_length_maxima(ngon_unit& PGs, const K_FLD::FloatArray& crd, E_Float& Lmin, E_Int& imin, E_Float& Lmax, E_Int& imax)
  {
    Lmin = K_CONST::E_MAX_FLOAT;
    Lmax = -1.;
    imin=imax=E_IDX_NONE;

    E_Int nb_pgs = PGs.size(), nb_nodes;
    E_Float L2;

    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      nb_nodes = PGs.stride(i);
      const E_Int* pNi = PGs.get_facets_ptr(i);

      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        L2 = K_FUNC::sqrDistance(crd.col(*(pNi + j) - 1), crd.col(*(pNi + (j + 1) % nb_nodes) - 1), 3);
        imin = (L2 < Lmin) ? i : imin;
        Lmin = (L2 < Lmin) ? L2 : Lmin;
        
        imax = (Lmax < L2) ? i : imax;
        Lmax = (Lmax < L2) ? L2 : Lmax;
      }
    }

    Lmin = ::sqrt(Lmin);
    Lmax = ::sqrt(Lmax);
  }

  struct link_t{
      link_t():count(-1){}
      E_Int nodes[2], count;
      bool add(E_Int n){ 
        if (count == -1)
        {
          nodes[++count] = n;
          return true;
        }
        else if (count == 0)
        {
          if (n != nodes[0])
            nodes[++count] = n;
          return true;
        }
        else // count = 1
        {
          if (n == nodes[0] || n == nodes[1])
            return true;
        }
        return false;//non manifold
      }
    };
  
  static E_Int flag_manifold_nodes(const ngon_t& ng, const K_FLD::FloatArray& crd, Vector_t<bool>& mnfld_nodes)
  {

    Vector_t<link_t> links(crd.cols());
    
    mnfld_nodes.resize(crd.cols(), true);
    E_Int nb_nodes;
    for (E_Int i = 0; i < ng.PGs.size(); ++i)
    {
      nb_nodes = ng.PGs.stride(i);
      const E_Int* ptr = ng.PGs.get_facets_ptr(i);

      for (E_Int j = 0; j < nb_nodes; ++j)
      {
        const E_Int Ni = *(ptr + j);
        const E_Int & Nip1 = *(ptr + (j + 1) % nb_nodes);

        mnfld_nodes[Ni - 1] = mnfld_nodes[Ni - 1] && links[Ni - 1].add(Nip1 - 1);
        mnfld_nodes[Nip1 - 1] = mnfld_nodes[Nip1 - 1] && links[Nip1 - 1].add(Ni - 1);
      }
    }

    return 0;

  }
  
  ///
  template <typename TriangulatorType>
  static E_Int __reorient_skin
(const K_FLD::FloatArray& coord, ngon_t& NGZ, const ngon_unit& pg_ext, const Vector_t<E_Int>& pg_ext_to_wPG, const ngon_unit& neighbors, Vector_t<E_Int>& orient)
{
  E_Int refPG(E_IDX_NONE), ori(E_IDX_NONE);
  E_Int err = __set_ref_PGs_for_orient<TriangulatorType>(coord, NGZ, refPG,ori);
      
  if (err)
    return err;
  
  if (refPG == E_IDX_NONE) // the mesh might be corrupted, anyway it must not happen
    return 1;
  //assert (refPG != E_IDX_NONE);

  assert (ori != E_IDX_NONE);
  
  //find the refPG
  E_Int iref = E_IDX_NONE;
  for (size_t i=0; i < pg_ext_to_wPG.size(); ++i)
  {
    if (pg_ext_to_wPG[i] == refPG)
    {
      iref = i;
      break;
    }
  }
  assert (iref != E_IDX_NONE);
  refPG=iref;
  
#ifdef DEBUG_NGON_T 
  /*K_FLD::IntArray tmpE;
  for (std::set<K_MESH::NO_Edge>::const_iterator it = walls.begin(); it != walls.end(); ++it)
    tmpE.pushBack(it->begin(), it->begin()+2);
  tmpE.shift(-1);
  std::cout << tmpE << std::endl;
  MIO::write("walls.mesh", coord, tmpE, "BAR");*/
#endif
    
  orient.resize(pg_ext.size(), 1);
  if (ori == -1)
    orient[iref]=-1;
  K_CONNECT::EltAlgo<K_MESH::Polygon>::reversi_connex(pg_ext, neighbors, refPG, orient);
  
  #ifdef DEBUG_NGON_T
  //Check orientation consistency
  {
    ngon_unit pgtmp = pg_ext;
    pgtmp.updateFacets();
  
    for (E_Int i = 0; i < pgtmp.size(); ++i)
    {
      if (orient[i] == -1)
      {
        E_Int s = pgtmp.stride(i);
        E_Int* p = pgtmp.get_facets_ptr(i);
        std::reverse(p, p + s);
      }
    }
    
    std::set<E_Int> faultys;
    E_Int er = K_CONNECT::EltAlgo<K_MESH::Polygon>::check_orientation_consistency(pgtmp, neighbors, faultys);
    if (faultys.empty())
    {
      Vector_t<E_Int> vfaultys, olids;
      vfaultys.insert(vfaultys.end(), faultys.begin(), faultys.end());
      NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs(coord, pgtmp, vfaultys);
      return er;
    }
  }
#endif
   
  return 0;
}
  
///
template <typename TriangulatorType>
static E_Int __set_ref_PGs_for_orient(const K_FLD::FloatArray& coord, ngon_t& ng, E_Int& PGref, E_Int& orient)
{
  E_Int err(0);
  K_FLD::ArrayAccessor<K_FLD::FloatArray> ac(coord);
  
  orient=E_IDX_NONE;
  PGref=E_IDX_NONE;
  
  for (E_Int PHi = 0; (PHi < ng.PHs.size()) && (PGref == E_IDX_NONE) && !err; ++PHi)
    if ((ng.PHs._type[PHi] == INITIAL_SKIN) || (ng.PHs._type[PHi] == CONNEXION_SKIN))
      err = __find_reference_orientation_triangle<TriangulatorType>(ng, coord, PHi, PGref, orient);
  
  return err;
}

///
template <typename TriangulatorType>
static E_Int __find_reference_orientation_triangle
(const ngon_t& ng, const K_FLD::FloatArray& coord, E_Int PHi, E_Int& PGgoal, E_Int& orient)
{
  PGgoal = E_IDX_NONE;
  orient = 0;
  
#ifdef DEBUG_NGON_T
  K_FLD::ArrayAccessor<K_FLD::FloatArray> ac(coord);
  NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::write("ng_for_ref_orient.plt", ac, ng);
  NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGT3s(coord, ng.PGs);
  //NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::highlight_PH<DELAUNAY::Triangulator>(ng, coord, PHi);
#endif
  
  // WARNING : assume here PHi is external
  E_Int nb_pg = ng.PHs.stride(PHi);
  const E_Int* ptr = ng.PHs.get_facets_ptr(PHi);
  
  Vector_t<E_Int> indices;
  indices.insert(indices.end(), ptr, ptr+nb_pg);
  for (size_t i = 0; i < indices.size(); ++i)indices[i] -= 1;//indcies start at one
  
  ngon_unit pgs;
  Vector_t<E_Int> oIds;
  ng.PGs.extract(indices, pgs, oIds);
  
  K_FLD::IntArray connectT3;
  Vector_t<E_Int> T3_to_nPG;
  E_Int err = ngon_t::triangulate_pgs<TriangulatorType>(pgs, coord, connectT3, T3_to_nPG);
  if (err)
    return err;
  
#ifdef DEBUG_NGON_T
  //MIO::write("PHref.mesh", coord, connectT3, "TRI", 0, &T3_to_nPG);
#endif
  
  {
  K_FLD::IntArray::const_iterator pS;
  E_Float Ni[3], P0[3], P1[3];
  E_Int dirp, dirm, PGi, PGj, PGprocessed;
  bool intersect;
  size_t sz = T3_to_nPG.size();
  for (size_t Ti = 0; Ti < sz; ++Ti)
  {
    // test only T3s belonging to skin PGs.
    PGi = oIds[T3_to_nPG[Ti]];
    if (!IS_EXTERNAL(ng.PGs._type[PGi]))
      continue;

 #ifdef DEBUG_NGON_T   
    /*std::vector<E_Int> colors;
    colors.resize(connectT3.cols(), 0);
    colors[Ti]=1;
    TRI_debug::write_wired("refT3wire.mesh", coord, connectT3, true, &colors);*/
#endif   
    
    pS = connectT3.col(Ti);
    K_MESH::Triangle::isoG(coord, pS, P0); // P0 is the center of Ti
    K_MESH::Triangle::normal(coord, pS, Ni); // Ni is the normal to Ti
    K_FUNC::sum<3>(P0, Ni, P1); //P0P1 == Ni
    
    dirp=dirm=0;
    PGprocessed = E_IDX_NONE;
    
    E_Float lambda, UV[2], min_d;
    E_Bool  coincident, parallel;
    
    for (size_t Tj = 0; Tj < sz; ++Tj)
    {
      PGj = oIds[T3_to_nPG[Tj]];
      if ((PGj == PGi) || (PGj ==  PGprocessed))
        continue;
      
      PGprocessed = E_IDX_NONE;
      const E_Int& N1 = connectT3(0,Tj);
      const E_Int& N2 = connectT3(1,Tj);
      const E_Int& N3 = connectT3(2,Tj);
      
      const E_Float* q0 = coord.col(N1);
      const E_Float* q1 = coord.col(N2);
      const E_Float* q2 = coord.col(N3);
      
      K_MESH::Triangle::planeLineMinDistance<3>(q0, q1, q2, P0, P1, E_EPSILON, true/*tol_is_absolute*/, lambda, UV, parallel, coincident, min_d);
      //intersect = K_MESH::Triangle::intersect<3>(q0, q1, q2, P0, P1, E_EPSILON, true, u0, u1, tx, overlap);
      //intersect |= (u0 != K_CONST::E_MAX_FLOAT) && (u1 == K_CONST::E_MAX_FLOAT);
      intersect= (((1. - UV[0] -  UV[1]) >= -E_EPSILON)  && (UV[0] >= -E_EPSILON) && (UV[1] >= -E_EPSILON));
      
      if (coincident)
        continue;
      else if (intersect)
      {
        if (lambda > E_EPSILON)++dirp;
        else if (lambda < -E_EPSILON)++dirm;
        PGprocessed = PGj;
      }
    }
    
    if (((dirp+dirm)%2 == 1))
    {
      if (dirp%2 == 0)orient = 1;
      else orient = -1;
      PGgoal = PGi;
      break;
    }
  }
  }
  return err;
}

///
template <typename TriangulatorType>
static E_Int reorient_skins(const TriangulatorType& t, const K_FLD::FloatArray& coord, ngon_t& wNG, bool& has_been_reversed)
{
  if (wNG.PHs.size()) assert(wNG.PHs.size() < wNG.PGs.size()); //i.e. not one ph per pg beacuse we need at least a closed PH

  has_been_reversed = false;

#ifdef FLAG_STEP
  chrono c;
  c.start();
#endif

  // 1. Get the skin (and connexion skin) PGs
  Vector_t<E_Int> oids, oids_cnx;
  ngon_unit pg_ext, pg_cnx;
  wNG.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);
  wNG.PGs.extract_of_type(CONNEXION_SKIN, pg_cnx, oids_cnx);

  oids.insert(oids.end(), oids_cnx.begin(), oids_cnx.end());
  pg_ext.append(pg_cnx);

  if (pg_ext.size() == 0)
    return 0;

#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__reorient_skins : extract_of_type : " << c.elapsed() << std::endl;
  c.start();
#endif

#ifdef DEBUG_NGON_T
  /*{
    //NGON_DBG::write("pgski1.plt", ACoordinate_t(coord), pg_ext);
    K_FLD::IntArray connectT3;
    Vector_t<E_Int> T3_to_nPG;
    E_Int err = ngon_t::triangulate_pgs<DELAUNAY::Triangulator>(pg_ext, coord, connectT3, T3_to_nPG);
    TRI_debug::write_wired("before_orient.mesh", coord, connectT3, true);
  }*/
#endif

  // 2. Build the neighbourhood for skin PGs
  ngon_unit neighbors;
  ngon_t::build_pg_neighborhood(pg_ext, neighbors);

  // update it to treat non-manifoldness as void neighbor.
  {
    // count edges sharing number
    std::map<K_MESH::NO_Edge, E_Int> edge_to_count;
    std::map<K_MESH::NO_Edge, E_Int>::iterator it;
    K_MESH::NO_Edge E;
    for (size_t i = 0; i < pg_ext.size(); ++i)
    {
      E_Int* pN = pg_ext.get_facets_ptr(i);
      E_Int nb_nodes = pg_ext.stride(i);

      for (size_t n = 0; n < nb_nodes; ++n)
      {
        E_Int ni = *(pN + n);
        E_Int nj = *(pN + (n + 1) % nb_nodes);
        E.setNodes(ni, nj);
        it = edge_to_count.find(E);
        if (it == edge_to_count.end())
          edge_to_count.insert(std::make_pair(E, 1));
        else
          ++it->second;
      }
    }
    // traverse the pgs and unset the non manifold neighbors
    for (size_t i = 0; i < pg_ext.size(); ++i)
    {
      E_Int* pN = pg_ext.get_facets_ptr(i);
      E_Int nb_nodes = pg_ext.stride(i);

      E_Int* pKn = neighbors.get_facets_ptr(i);

      for (size_t n = 0; n < nb_nodes; ++n)
      {
        E_Int ni = *(pN + n);
        E_Int nj = *(pN + (n + 1) % nb_nodes);
        E.setNodes(ni, nj);

        if (edge_to_count[E] != 2)
          *(pKn + n) = E_IDX_NONE;
      }
    }
  }

#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__reorient_skins : build_pg_neighborhood : " << c.elapsed() << std::endl;
  c.start();
#endif

  // Color to get connex parts
  E_Int nb_connex = 1;
  Vector_t<E_Int> colors;
  if (wNG.PHs.size() > 1)
  {
    K_CONNECT::EltAlgo<K_MESH::Polygon>::coloring(neighbors, colors);
    nb_connex = 1 + *std::max_element(colors.begin(), colors.end());
  }

#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__reorient_skins : coloring : " << c.elapsed() << std::endl;
  c.start();
#endif

  E_Int err = 0;
  Vector_t<E_Int> orient;
  if (nb_connex > 1)
  {
    Vector_t<bool> keep;
    ngon_t NGZ;
    for (E_Int z = 0; z < nb_connex; ++z)
    {
#ifdef FLAG_STEP
      c.start();
#endif

      NGZ = wNG;
      NGZ.PHs._type = wNG.PHs._type;
      NGZ.PGs._type = wNG.PGs._type;

      keep.resize(NGZ.PGs.size(), false);
      for (size_t i = 0; i < colors.size(); ++i) keep[oids[i]] = (colors[i] == z);
      NGZ.keep_PHs_having_PGs(keep);

#ifdef FLAG_STEP
      if (chrono::verbose >1) std::cout << "__reorient_skins : prepare connex part number : " << z << " " << c.elapsed() << std::endl;
      c.start();
#endif

      err = __reorient_skin<TriangulatorType>(coord, NGZ, pg_ext, oids, neighbors, orient); //reorient a connex part NGZ of wNG. pg_ext are related to wNG
      if (err)
      {
        std::cout << "could not orient properly zone " << z << std::endl;
        
#ifdef DEBUG_NGON_T
        std::ostringstream o;
        o << "reorient_error_zone_" << z << ".tp";
        K_FLD::IntArray cnto;
        NGZ.export_to_array(cnto);
        MIO::write(o.str().c_str(), coord, cnto, "NGON");
#endif
        return err;
      }

#ifdef FLAG_STEP
      if (chrono::verbose >1) std::cout << "__reorient_skins : __reorient_skin number : " << z << " " << c.elapsed() << std::endl;
      c.start();
#endif
    }
  }
  else //one connex zone
  {
    err = __reorient_skin<TriangulatorType>(coord, wNG, pg_ext, oids, neighbors, orient); //reorient a connex part NGZ of wNG. pg_ext are related to wNG

#ifdef FLAG_STEP
    if (chrono::verbose >1) std::cout << "__reorient_skins : __reorient_skin : " << c.elapsed() << std::endl;
    c.start();
#endif
  }

  //Apply new orientation
  for (E_Int i = 0; i < pg_ext.size(); ++i)
  {
    if (orient[i] == -1)
    {
      E_Int PGi = oids[i];
      E_Int s = wNG.PGs.stride(PGi);
      E_Int* p = wNG.PGs.get_facets_ptr(PGi);
      std::reverse(p, p + s);
      has_been_reversed = true;
    }
  }

#ifdef FLAG_STEP
  if (chrono::verbose >1) std::cout << "__reorient_skins : apply orientation  : " << c.elapsed() << std::endl;
#endif

#ifdef DEBUG_NGON_T
  /*{
    K_FLD::IntArray connectT3;
    Vector_t<E_Int> T3_to_nPG;
        
    Vector_t<E_Int> oids, oids_cnx;
    ngon_unit pg_ext, pg_cnx;
    wNG.PGs.extract_of_type(INITIAL_SKIN, pg_ext, oids);
    wNG.PGs.extract_of_type(CONNEXION_SKIN, pg_cnx, oids_cnx);

    oids.insert(oids.end(), oids_cnx.begin(), oids_cnx.end());
    pg_ext.append(pg_cnx);
    
    E_Int err = ngon_t::triangulate_pgs<DELAUNAY::Triangulator>(pg_ext, coord, connectT3, T3_to_nPG);
    TRI_debug::write_wired("oriented.mesh", coord, connectT3, true);
  }*/
#endif

  return 0;
}

///
static E_Int reorient_connex_PGs(ngon_unit& PGs, bool reverse_first)
{
  // WARNING : use the first Polygon as the reference. Caller can reverse it if required.
  
  PGs.updateFacets();
  
  Vector_t<E_Int> orient(PGs.size(), 1);
  if (reverse_first)
    orient[0]=-1;
  
  ngon_unit neighbors;
  build_pg_neighborhood(PGs, neighbors);
  
  K_CONNECT::EltAlgo<K_MESH::Polygon>::reversi_connex(PGs, neighbors, 0/*reference PG*/, orient);
  
  //Apply new orientation
  for (E_Int PGi = 0; PGi < PGs.size(); ++PGi)
  {
    if (orient[PGi] == -1)
    {
      E_Int s = PGs.stride(PGi);
      E_Int* p = PGs.get_facets_ptr(PGi);
      std::reverse(p, p + s);
    }
  }
  
  return 0;
}
  
/// coloring-like algo
void
build_F2E(const ngon_unit& neighbors, K_FLD::IntArray& F2E) const
{

  // WARNING : assuming that external PGs have been propoerly oriented toward exterior
  E_Int NB_PHS(PHs.size());
  E_Int NB_PGS(PGs.size());

  F2E.clear();
  F2E.resize(2, NB_PGS, E_IDX_NONE);

  Vector_t<E_Int> exPH(NB_PHS, E_IDX_NONE);

  //init for external elements : setting the left
  E_Int nb_phs = neighbors.size();
  for (size_t PHi = 0; PHi < nb_phs; ++PHi)
  {
    const E_Int* pGi = PHs.get_facets_ptr(PHi);
    const E_Int& nb_neigh = PHs.stride(PHi);
    const E_Int* pKn = neighbors.get_facets_ptr(PHi);

    for (size_t i = 0; (i<nb_neigh); ++i)
      if (*(pKn + i) == E_IDX_NONE){
        //std::cout << "PHi/PGi : " << PHi << "/" << *(pGi + i) - 1 << std::endl;
        F2E(0, *(pGi + i) - 1) = PHi;  // left PH
        exPH[PHi] = *(pGi + i); //store i-th external PG for this external PH
        break;
      }
  }

  E_Int PHseed(0), PH, iref, PGref;
  Vector_t<bool> processed(NB_PHS, false);
  ngon_unit PGneighbors, pgs;
  Vector_t<E_Int> oids, orient, cpool;
  bool revers;

  while (1) // as many iter as connex bits
  {
    while ((PHseed < NB_PHS) && ((processed[PHseed] == true) || (exPH[PHseed] == E_IDX_NONE))) ++PHseed; // take first external element
    if (NB_PHS - 1 < PHseed)
      return;

    cpool.push_back(PHseed);

    while (!cpool.empty())
    {
      PH = cpool.back();
      cpool.pop_back();

      if (processed[PH])
        continue;

      processed[PH] = true;

      const E_Int* pGi = PHs.get_facets_ptr(PH);
      const E_Int& nb_neigh = PHs.stride(PH);
      const E_Int* pKn = neighbors.get_facets_ptr(PH);

      // build PG neighborhood for that given PH
      PGs.extract(pGi, nb_neigh, pgs, oids);
      build_pg_neighborhood(pgs, PGneighbors);

      orient.clear();
      orient.resize(nb_neigh, 1);
      revers = false;
      PGref = exPH[PH]; //satrt at 1, can be negative
      if (PGref < 0)
      {
        revers = true;
        PGref = -PGref;
      }

      //found it in oids
      iref = E_IDX_NONE;
      for (size_t i = 0; i < oids.size(); ++i)
        if (PGref == (oids[i] + 1)){ iref = i; break; }
      assert(iref != E_IDX_NONE);
      if (revers)
        orient[iref] = -1;

      K_CONNECT::EltAlgo<K_MESH::Polygon>::reversi_connex(pgs, PGneighbors, iref, orient);

      for (size_t i = 0; (i<nb_neigh); ++i)
      {
        E_Int PGi = *(pGi + i);
        const E_Int& PHn = *(pKn + i);

        F2E(0, PGi - 1) = PH;
        F2E(1, PGi - 1) = PHn;

        if (PHn == E_IDX_NONE)
          continue;

        exPH[PHn] = -PGi; //starts at 1 so logic with negative value is OK

        if (orient[i] == -1){ //swap
          std::swap(F2E(0, PGi - 1), F2E(1, PGi - 1));
          exPH[PHn] = PGi;
        }

        if (!processed[PHn])
          cpool.push_back(PHn);
      }
    }
  }
}

#ifdef DEBUG_BOOLEAN
// Writes a single zone into a file (element type must be specified to distinguish QUAD/TETRA, it is not ambiguous otherwise.)
static E_Int write
(const char* filename, const K_FLD::FloatArray& crd, const ngon_t& ng, const std::vector<bool>* mask=0, const std::vector<E_Int>* colors=0)
{
  K_FLD::IntArray cnt;
  ng.export_to_array(cnt);
  return MIO::write(filename, crd, cnt, "NGON", mask, colors);
}
#endif

///
/*E_Int keep_PHs(const Vector_t<E_Int>& PHs_tokeep)
{
  ngon_unit keptPHs;
  for (size_t i = 0; i < PHs_tokeep.size(); ++i)
    keptPHs.__add(PHs, PHs_tokeep[i]);

  PHs=keptPHs;
  PHs.updateFacets();

  return 0;
}
  
///
E_Int remove_PHs(const Vector_t<E_Int>& PHs_toremove)
{
  std::set<E_Int> toremset(PHs_toremove.begin(), PHs_toremove.end());
  ngon_unit keptPHs;
 
  E_Int sz = PHs.size();
  for (size_t i = 0; i < sz; ++i)
  {
    if (toremset.empty() || (toremset.find(i) == toremset.end()))
      keptPHs.__add(PHs, i);
    else
      toremset.erase(i);
  }
      
  PHs=keptPHs;
  PHs.updateFacets();
    
  return 0;
}
  
///
E_Int reorient_PGs(const Vector_t<E_Int>& PGs_toreverse)
{
  E_Int sz1(PGs_toreverse.size()), sz2;
  for (size_t i = 0; i < sz1; ++i)
  {
    const E_Int& PGi = PGs_toreverse[i];
    sz2 = PGs.get_facets_nb(PGi) >> 1; //half size
    //do the reverse by swapping
    E_Int* ptr = PGs.get_facets_ptr(PGi);
    for (size_t j = 0; j < sz2; ++j)
      std::swap(ptr[j], ptr[sz2-1-j]);
  }
  return 0;
}*/
    
  ngon_unit PGs;
  ngon_unit PHs;
  
};

#endif	/* NGON_T_HXX */

