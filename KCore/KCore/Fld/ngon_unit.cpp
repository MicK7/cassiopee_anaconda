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
#include "ngon_unit.h"

#include "Connect/connect.h"
#include "Connect/IdTool.h"
#include <assert.h>
#include <algorithm>
#include <map>

/// Methods Definitions
ngon_unit::ngon_unit(const E_Int* begin):_dirty(true)
{
  E_Int l = *(begin+1)+2;
  const E_Int* end = begin+l;
  _NGON.clear();
  _NGON.insert(_NGON.end(), begin, end);
}

///
void ngon_unit::updateFacets() const
{
  if (_dirty && !_NGON.empty())
  {
    K_CONNECT::getPosFacets(&_NGON[0], 2, _NGON[0], _facet);
    _dirty=false;
  }
}

///
ngon_unit& ngon_unit::operator=(const ngon_unit& ng)
{
  _NGON.clear();
  _NGON = ng._NGON;
  _type = ng._type;
  _ancEs = ng._ancEs;
  _dirty=ng._dirty;
  if (!_dirty)
    _facet = ng._facet;
  return *this;
}

///
ngon_unit::ngon_unit(const ngon_unit& ngin)
{
  *this = ngin;
}

///
bool ngon_unit::is_consistent()
{
  if (!_type.empty() && _type.size() != _facet.size())
    return false;
  if (_ancEs.cols() != 0 && _ancEs.cols() != _facet.size())
    return false;
  return true;
}


///
void ngon_unit::add(E_Int n, const E_Int* facet_ptr)
{
  // molecule here is one PH or PG
  if (_NGON.empty())
  {
    _NGON.resize(2, 0);
    _dirty = false;
  }

  _NGON[0]+=1;
  _NGON.push_back(n);
  _NGON.insert(_NGON.end(), facet_ptr, facet_ptr + n);
  _NGON[1] = _NGON.size() - 2;

  if (!_dirty)
    _facet.push_back(_NGON.size() - n - 1);
}

///
E_Int ngon_unit::get_facets_max_id() const
{
  updateFacets();
  
  E_Int nb_nodes, nb_pg(size()), id;  
  E_Int maxId=-1;
  // Loop through the elements and increment face_count
  for (E_Int i = 0; i < nb_pg; ++i)
  {
    nb_nodes = stride(i);
    for (E_Int j = 0; j < nb_nodes; ++j)
    {
      id = get_facet(i,j);
      maxId = (maxId < id) ? id : maxId;
    }
  }
  return maxId;
}

///
void ngon_unit::append(const ngon_unit& cngon_in)
{
  if (cngon_in.size() == 0)
    return;

  if (_NGON.empty())
    *this=cngon_in;
  else
  {
    E_Int sz1 = _NGON.size();
    _NGON[0] += cngon_in._NGON[0];
    _NGON.insert(_NGON.end(), cngon_in._NGON.begin()+2, cngon_in._NGON.end());
    _NGON[1]=_NGON.size()-2;
    if (!_type.empty()) //transfer the info if it is already computed for the left hand side
      _type.insert(_type.end(), cngon_in._type.begin(), cngon_in._type.end());
    if (_ancEs.cols()) //transfer the info if it is already computed for the left hand side
      _ancEs.pushBack(cngon_in._ancEs);
      
    if (!_dirty && !cngon_in._dirty)
    {
      E_Int sz2 = _facet.size();
      _facet.insert(_facet.end(), cngon_in._facet.begin(), cngon_in._facet.end());
      K_CONNECT::IdTool::shift(_facet, sz2, sz1-2); //tail shift
    }
    else
      _dirty=true;
  }
}

///
void ngon_unit::append(const Vector_t<E_Int>& NGON)
{
  if (_NGON.empty())
    _NGON=NGON;
  else
  {
    _NGON[0] += NGON[0];
    _NGON.insert(_NGON.end(), NGON.begin()+2, NGON.end());
    _NGON[1]=_NGON.size()-2;
  }
  _dirty=true;
}

///
void ngon_unit::__add (const ngon_unit& ng, E_Int ith)
{
  if (_NGON.empty())
  {
    _NGON.resize(2,0);
    _dirty=false;
  }
  
  _NGON[0]++;
  E_Int sz = ng.stride(ith);
  const E_Int* ptr = ng.get_facets_ptr(ith);
  _NGON.push_back(sz);
  _NGON.insert(_NGON.end(), ptr, ptr+sz);
  _NGON[1]=_NGON.size()-2;
  if (!_dirty )
    _facet.push_back(_NGON.size()-1-sz);
  if (!ng._type.empty()) //transfer externality if the info is already computed for the left hand side. fixme hpc
    _type.push_back(ng._type[ith]);
  if (ng._ancEs.cols()) //transfer ancestors if the info is already computed for the left hand side. fixme hpc
    _ancEs.pushBack(ng._ancEs.col(ith), ng._ancEs.col(ith)+2);
}

///
E_Int ngon_unit::remove_entities (const Vector_t<E_Int>& to_remove, Vector_t<E_Int>& nids) //0-based nids
{
  nids.clear();

  if (to_remove.empty())
    return 0;

  E_Int sz(size()), count(0);
  std::set<E_Int> remset(to_remove.begin(), to_remove.end());//fixme
  std::set<E_Int>::const_iterator remend(remset.end());

  nids.resize(sz, E_IDX_NONE);

  updateFacets();

  ngon_unit ngu;

  for (E_Int i = 0; (i < sz); ++i)
  {
    if (remset.find(i) != remend)
      continue;

    ngu.__add(*this, i);
    nids[i]=count++;
  }
  
  if (!_type.empty())
  {
    ngu._type=_type;
    if (count != sz) //nothing has been removed
      K_CONNECT::IdTool::compact(ngu._type, nids);
  }
  if (_ancEs.cols())
  {
    ngu._ancEs=_ancEs;
    if (count != sz) //nothing has been removed
      K_CONNECT::IdTool::compact(ngu._ancEs, nids);
  }

  *this=ngu;
  updateFacets();
  
#ifdef DEBUG_NGON_UNIT
  assert(this->is_consistent());
#endif
  
  return sz-size();
}

///
void ngon_unit::remove_facets(const Vector_t<E_Int>& facet_ids, Vector_t<E_Int>& nids) //0-based nids
{

  ngon_unit ngu;
  Vector_t<E_Int> molecule;
  
  updateFacets();
  
  nids.clear();
  nids.resize(size(), E_IDX_NONE);

  E_Int sz, count(-1);
  //std::cout << "nb elts : " << PHs.size() << std::endl;
  for (E_Int i = 0; i < size(); ++i)
  {
    molecule.clear();
    sz = stride(i);
    if (sz == 0) //fixme : required ?
      continue;
    molecule.push_back(sz);
    for (E_Int j = 0; j < sz; ++j)
    {
      const E_Int& Fi = get_facet(i, j);
      /*if (Fi < 0 || Fi > facet_ids.size()-1)
      {
      std::cout << "elt : " <<  i << std::endl;
      std::cout << "elt stride : " << sz << std::endl;
      std::cout << "BOUM : " << Fi << std::endl;
      return 1;
      }*/
      if (facet_ids[Fi-1] == E_IDX_NONE)
        continue;
      molecule.push_back(facet_ids[Fi-1]+1);
    }
    molecule[0] = molecule.size() - 1;
    if (molecule[0]>0)//i.e has at least one valid face.
    {
      ngu.add(molecule);
      nids[i]=++count;
    }
  }
  
  //transfer externality
  ngu.__transfer_attributes(*this, nids);

  *this = ngu;
}
  
///
void ngon_unit::shift(E_Int val)
{
  if (val == 0) return;

  updateFacets();

  size_t sz(_facet.size()), nb_nodes;  
  for (size_t i = 0; i < sz; ++i)
  {
    nb_nodes = stride(i);
    E_Int* ptr = get_facets_ptr(i);
    for (size_t j = 0; j < nb_nodes; ++j)
      *(ptr+j) += val;
  }
}

///set all facet values to E_IDX_NONE (useful for building a neighbor table)
void
ngon_unit::reset_facets()
{
  updateFacets();
  
  size_t sz(_facet.size()), nb_nodes;  
  for (size_t i = 0; i < sz; ++i)
  {
    nb_nodes = stride(i);
    E_Int* ptr = get_facets_ptr(i);
    for (size_t j = 0; j < nb_nodes; ++j)
      *(ptr+j) = E_IDX_NONE;
  }
}

///
void ngon_unit::extract
(const Vector_t<E_Int>& indices, ngon_unit& ng_out, Vector_t<E_Int>& oldIds) const
{
  // 0-based indices
  ng_out.clear();
  oldIds.clear();
  oldIds.reserve(indices.size());
  updateFacets();
  
  for (size_t i = 0; i < indices.size(); ++i)
  {
    ng_out.__add(*this, indices[i]);
    oldIds.push_back(indices[i]);
  }
}

///
void ngon_unit::extract
(const E_Int* ptr, E_Int n, ngon_unit& ng_out, Vector_t<E_Int>& oids) const
{
  // 0-based indices
  ng_out.clear();
  oids.clear();
  oids.reserve(n);
  updateFacets();
  
  for (E_Int i = 0; i < n; ++i)
  {
    ng_out.__add(*this, *(ptr+i)-1);
    oids.push_back(*(ptr+i)-1);
  }
}

///
void ngon_unit::extract_of_type (E_Int FLAG, ngon_unit& ng_out, Vector_t<E_Int>& oldIds) const
{
  Vector_t<E_Int> indices;
  get_indices_of_type(FLAG, indices);
  ng_out.clear();
  extract(indices, ng_out, oldIds);//get E32 skin faces
}

///
void ngon_unit::get_indices_of_type(E_Int FLAG, Vector_t<E_Int>& indices) const
{
  // 0-based indices
  for (size_t i = 0; i < _type.size(); ++i)
    if (_type[i] == FLAG)indices.push_back(i);
}

///
void ngon_unit::flag_indices_of_type (E_Int FLAG, Vector_t<bool>& flag) const
{

  flag.clear();
  flag.resize(_type.size(), false);
  for (size_t i = 0; i < _type.size(); ++i) flag[i] = (_type[i] == FLAG);
}

///
void ngon_unit::find_elts_with_facets(const std::set<E_Int>& facets, std::vector<E_Int> & elts) const
{
  elts.clear();
  
  updateFacets();
  
  E_Int nb_facets, nb_elts(size()), id;
  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      if (facets.find(id) != facets.end())
      {
        elts.push_back(i);
        break;
      }
    }
  }  
}

///
void ngon_unit::unique_indices(std::vector<E_Int>& indices) const
{
  std::set<E_Int> tmp;
  E_Int nb_facets, nb_elts(size()), id;
  
  indices.clear();
  
  updateFacets();
  
  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      tmp.insert(id);
    }
  }
  
  indices.insert(indices.end(), tmp.begin(), tmp.end());
}

/// NOT CONSECUTIVE : (N0, N1, N1, N2, N3, N3, N3, N4, N1, N7...) => (N0, N1, N2, N3, N4, N7...)
//any doublon occuring is removed => not ok when facets are node (PG destructive), ok when facet is PG
E_Int ngon_unit::remove_duplicated()
{
  std::set<E_Int> stmp;
  E_Int nb_facets, nbf, nb_elts(size()), id, nbe(0);
  ngon_unit ngtmp;
  
  updateFacets();
  
  ngtmp.clear(); 
  ngtmp._NGON.resize(2, 0);

  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    
    if (nb_facets == 1)
    {
      ngtmp._NGON.push_back(1);
      ngtmp._NGON.push_back(get_facet(i,0));
      ++nbe;
      continue;
    }
    
    stmp.clear();
    nbf=0;
    ngtmp._NGON.push_back(0);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      //std::cout << "i / id " << i << " : " << id << std::endl;
      if (stmp.insert(id).second)
      {
        ngtmp._NGON.push_back(id);
        ++nbf;
      } 
    }
    
    if (nbf)
    {
      ++nbe;
      ngtmp._NGON[ngtmp._NGON.size()-1-nbf]=nbf;//set the nb of facets
    }
    else
      ngtmp._NGON.pop_back();
  }
  
  ngtmp._NGON[0]=nbe;
  ngtmp._NGON[1]=ngtmp._NGON.size()-2;
  
  E_Int prev_sz = this->_NGON[1];
  
  // hack to preserve info (type and history) : must work because this function does not change the number of facets
  if (!_type.empty())
    ngtmp._type=this->_type;
  if (_ancEs.cols())
    ngtmp._ancEs=this->_ancEs;
  
  *this=ngtmp;
  
  return this->_NGON[1] - prev_sz;
}

/// CONSECUTIVE (N0, N1, N1, N2, N3, N3, N3, N4, N1, N7...) => (N0, N1, N2, N3, N4, N1, N7...) :
//=> ok when facets are node (PG safe), not ok when facet is PG because some duplicates remain
E_Int ngon_unit::remove_consecutive_duplicated()
{
  E_Int nb_facets, nbf, nb_elts(size()), id, prev_id, nbe(0);
  ngon_unit ngtmp;
  
  updateFacets();
  
  ngtmp.clear(); 
  ngtmp._NGON.resize(2, 0);
  
  std::vector<E_Int> nids(nb_elts, E_IDX_NONE);
  
  // 
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    
    if (nb_facets == 1)
    {
      ngtmp._NGON.push_back(1);
      ngtmp._NGON.push_back(get_facet(i,0));
      nids[i] = nbe++;
      continue;
    }
    
    nbf=0;
    ngtmp._NGON.push_back(0);
    prev_id=get_facet(i,nb_facets-1);//initialized with the last one to see if not equal to the first one
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      id = get_facet(i,j);
      //std::cout << "i / id " << i << " : " << id << std::endl;
      if (id != prev_id)
      {
        ngtmp._NGON.push_back(id);
        prev_id=id;
        ++nbf;
      } 
    }
    if (nbf)
    {
      nids[i] = nbe++;
      ngtmp._NGON[ngtmp._NGON.size()-1-nbf]=nbf;//set the nb of facets
    }
    else
      ngtmp._NGON.pop_back();
  }
  
  ngtmp._NGON[0]=nbe;
  ngtmp._NGON[1]=ngtmp._NGON.size()-2;
  
  E_Int prev_sz = (this->_NGON[1]);
  
  //transfer attributes
  ngtmp.__transfer_attributes(*this, nids);
  
  *this=ngtmp;
  
  return ngtmp._NGON[1] - prev_sz;
}

///
/*
Une face est degeneree quand le nbre de pts identiques est superieur
 au nombre de noeuds-2.
 nvert=1 -> 1D -> jamais degenere
 nvert=2 -> 2D -> egalite=1 (degenere)
 nvert=3 -> 3D (face TRI) -> egalite >= 1 (degenere)
 nvert=4 -> 3D (face QUAD) -> egalite >= 2 (degenere)
 nvert=5 -> 3D -> egalite >= 3 (degenere)
*/
//=============================================================================
void ngon_unit::get_degenerated(Vector_t<E_Int>& indices) 
{
  updateFacets();
  
  E_Int s, nb_elts(size()), eq;
  std::set<E_Int> unic;
  indices.clear();
  
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    s = stride(i);
    if (s==1) continue;
    unic.clear();
    unic.insert(get_facets_ptr(i), get_facets_ptr(i)+s);
    eq=unic.size();//s - unic.size()
    eq = (s==2) ? eq+1 : eq;
    if (eq <= 2) // i.e s-2 <= s-u
      indices.push_back(i); 
  }
}

///
void ngon_unit::change_indices (const Vector_t<E_Int>& nIds, bool zerobased)
{
  updateFacets();
  
  E_Int nb_facets, nb_elts(size());
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    nb_facets = stride(i);
    for (E_Int j = 0; j < nb_facets; ++j)
    {
      E_Int& id = get_facet(i,j);
      if (zerobased)
        id = nIds[id-1]+1;
      else
        id = nIds[id];
    }
  }
}

bool ngon_unit::is_fixed_stride(E_Int& strd)
{
  updateFacets();
  
  strd=stride(0);
  
  E_Int nb_elts(size());
  for (E_Int i = 1; i < nb_elts; ++i)
  {
    if (strd != stride(i))
      return false;
  }
  return true;
}

/// convert a fixed-stride storage to a ngon_unit storage
void ngon_unit::convert_fixed_stride_to_ngon_unit(const K_FLD::IntArray&cnt, ngon_unit& nguo)
{
  //
  E_Int nb_elts = cnt.cols();
  E_Int stride = cnt.rows();

  E_Int sz = (stride + 1)*nb_elts + 2;
  nguo._NGON.resize(sz);

  nguo._NGON[0] = nb_elts;
  nguo._NGON[1] = sz - 2;

  E_Int pos = 2;
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    const E_Int* p = cnt.col(i);
    nguo._NGON[pos++] = stride;
    for (E_Int j = 0; j < stride; ++j)
      nguo._NGON[pos++] = *(p + j);
  }

  nguo.updateFacets();
}

void ngon_unit::create_pgs_changes_graph(const Vector_t<E_Int>& oids, bool zero_based_in, E_Int appendshift, ngon_unit& nguo) const 
{
  std::map<E_Int, std::vector<E_Int> > molecules;

  E_Int sz = oids.size();
  E_Int nb_pgs = this->size();

  E_Int shift = (zero_based_in) ? 0 : -1;
  E_Int max_id = -1;

  for (size_t i = 0; i < oids.size(); ++i)
  {
    E_Int id = oids[i] + shift;
    molecules[id].push_back(i + 1 + appendshift);
    max_id = (max_id < id) ? id : max_id;
  }

  assert(max_id < nb_pgs);

  E_Int empty = E_IDX_NONE;
  Vector_t<E_Int> molecule;
  std::map<E_Int, std::vector<E_Int> >::const_iterator it;

  nguo.clear();
  nguo._NGON.resize(2, 0);

  for (E_Int i = 0; i < nb_pgs; ++i)
  {
   /* if (toremove[i])
    {
      nguo._NGON.push_back(0);
      ++nguo._NGON[0];
      nguo._NGON[1] = nguo._NGON.size() - 2;
      continue;
    }*/
    it = molecules.find(i);
    if (it == molecules.end())
      nguo.add(1, &empty);
    else
      nguo.add(it->second.size(), &it->second[0]);
  }

  nguo.updateFacets();
}

/// transfer attribute
void ngon_unit::__transfer_attributes(const ngon_unit& ngi, const Vector_t<E_Int>& nids)
{
  if (ngi._type.empty() && (ngi._ancEs.cols() == 0)) return;
  
  //transfer externality
  if (!ngi._type.empty())
  {
    _type = ngi._type;
    K_CONNECT::IdTool::compact(_type, nids);
  }
  if (ngi._ancEs.cols() != 0)
  {
    _ancEs = ngi._ancEs;
    K_CONNECT::IdTool::compact(_ancEs, nids);
  }
  
  this->updateFacets();

//#ifdef DEBUG_NGON_UNIT
    assert(is_consistent());
//#endif
}

///
void ngon_unit::__transfer_attributes(const ngon_unit& ngi, const Vector_t<bool>& keep)
{
  if (ngi._type.empty() && (ngi._ancEs.cols() == 0)) return;
  
  //transfer externality
  K_CONNECT::keep pred(keep);
  
  if (!ngi._type.empty())
  {
    _type = ngi._type;
    K_CONNECT::IdTool::compress(_type, pred);
  }
  if (ngi._ancEs.cols() != 0)
  {
    _ancEs = ngi._ancEs;
    K_CONNECT::IdTool::compress(_ancEs, pred);
  }
  
  this->updateFacets();
  
//#ifdef DEBUG_NGON_UNIT
    assert(is_consistent());
//#endif
}
