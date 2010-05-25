/**********************************************************************
generic.h - Handle generic data classes. Custom data for atoms, bonds, etc.
 
Copyright (C) 1998-2001 by OpenEye Scientific Software, Inc.
Some portions Copyright (C) 2001-2006 by Geoffrey R. Hutchison
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#ifndef OB_GENERIC_H
#define OB_GENERIC_H

#include <openbabel/babelconfig.h>

#include <string>
#include <vector>
#include <map>

#include <openbabel/math/spacegroup.h>
#include <openbabel/obutil.h>
#include <openbabel/base.h>

namespace OpenBabel
{

  // Forward declarations
  class OBBase;
  class OBAtom;
  class OBBond;
  class OBMol;
  class OBRing;

  //! \class OBCommentData generic.h <openbabel/generic.h>
  //! \brief Used to store a comment string (can be multiple lines long)
 class OBAPI OBCommentData : public OBGenericData
  {
  protected:
    std::string _data;
  public:
    OBCommentData();
    OBCommentData(const OBCommentData&);
    virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new OBCommentData(*this);}
		
    OBCommentData& operator=(const OBCommentData &src);

    void          SetData(const std::string &data)
    { _data = data; Trim(_data); }
    void          SetData(const char *d)
    {_data = d; Trim(_data);     }
    const std::string &GetData()              const
    {        return(_data);      }
    virtual const std::string &GetValue()              const  
    {        return(_data);      }
  };

  //! \class OBExternalBond generic.h <openbabel/generic.h>
  //! \brief Used to store information on an external bond 
  //! (e.g., SMILES fragments)
  class OBAPI OBExternalBond
  {
    int     _idx;
    OBAtom *_atom;
    OBBond *_bond;
  public:
  OBExternalBond(): _idx(0), _atom(NULL), _bond(NULL) {}
    OBExternalBond(OBAtom *,OBBond *,int);
    OBExternalBond(const OBExternalBond &);
    ~OBExternalBond()   {}

    int     GetIdx()  const    {        return(_idx);    }
    OBAtom *GetAtom() const    {        return(_atom);   }
    OBBond *GetBond() const    {        return(_bond);   }
    void SetIdx(int idx)       {        _idx = idx;      }
    void SetAtom(OBAtom *atom) {        _atom = atom;    }
    void SetBond(OBBond *bond) {        _bond = bond;    }
  };

  //! \class OBExternalBondData generic.h <openbabel/generic.h>
  //! \brief Used to store information on external bonds (e.g., in SMILES fragments)
 class OBAPI OBExternalBondData : public OBGenericData
  {
  protected:
    std::vector<OBExternalBond> _vexbnd;
  public:
    OBExternalBondData();
		
    //Copying is not used and too much work to set up
    virtual OBGenericData* Clone(OBBase* /*parent*/) const{return NULL;}
    
    void SetData(OBAtom*,OBBond*,int);
    std::vector<OBExternalBond> *GetData()
      {
        return(&_vexbnd);
      }
  };

  //! \class OBPairData generic.h <openbabel/generic.h>
  //! \brief Used to store arbitrary text attribute/value relationships.
  //!
  //! Ideal for arbitrary text descriptors for molecules, atoms, bonds, residues,
  //!  e.g. in QSAR.
 class OBAPI OBPairData : public OBGenericData
  {
  protected:
    std::string _value; //!< The data for this key/value pair
  public:
    OBPairData();
    virtual OBGenericData* Clone(OBBase* /*parent*/) const
      {return new OBPairData(*this);}
    void    SetValue(const char *v)        {      _value = v;    }
    void    SetValue(const std::string &v) {      _value = v;    }
    virtual const std::string &GetValue() const
    {      return(_value);    }
  };

  //! \class OBPairTemplate generic.h <openbabel/generic.h>
  //! \brief Used to store arbitrary attribute/value relationsips of any type.
  // More detailed description in generic.cpp
  template <class ValueT>
    class OBAPI OBPairTemplate : public OBGenericData
  {
  protected:
    ValueT _value; //!< The data for this key/value pair
  public:
  OBPairTemplate():
    OBGenericData("PairData", OBGenericDataType::PairData) {};
    void SetValue(const ValueT t)             { _value = t;     }
    virtual const ValueT &GetGenericValue() const    { return(_value); }
  };

  //! Store arbitrary key/value integer data like OBPairData
  typedef OBPairTemplate<int>     OBPairInteger;
  //! Store arbitrary key/value floating point data like OBPairData
  typedef OBPairTemplate<double>  OBPairFloatingPoint;

  //! \class OBSetData generic.h <openbabel/generic.h>
  //! \brief Used to store arbitrary attribute/set relationships.
  //! Should be used to store a set of OBGenericData based on an attribute.
 class OBAPI OBSetData : public OBGenericData
  {
  protected:
    std::vector<OBGenericData *> _vdata;
  public:
  OBSetData() : OBGenericData("SetData", OBGenericDataType::SetData) {}
    virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new OBSetData(*this);}

    //! Add an OBGenericData element to the set.
    void AddData(OBGenericData *d)
    {
      if(d)
        {
          _vdata.push_back(d);
        }
    }

    //! Set the array of data to a new vector
    void SetData(std::vector<OBGenericData *> &vdata)
    {
      _vdata = vdata;
    }

    //! \return the OBGenericData associate with the attribute name parameter.
    OBGenericData *GetData(const char *s)
    {
      std::vector<OBGenericData*>::iterator i;

      for (i = _vdata.begin();i != _vdata.end();++i)
        if ((*i)->GetAttribute() == s)
          return(*i);

      return(NULL);
    }

    //! \return the OBGenericData associate with the attribute name parameter.
    OBGenericData *GetData(const std::string &s)
    {
      std::vector<OBGenericData*>::iterator i;

      for (i = _vdata.begin();i != _vdata.end();++i)
        if ((*i)->GetAttribute() == s)
          return(*i);

      return(NULL);
    }

    //! Gets the entire set.
    virtual const std::vector<OBGenericData *> &GetData() const //now virtual and const
    {
      return(_vdata);
    }

    //! Get the begin iterator.
    std::vector<OBGenericData*>::iterator GetBegin()
      {
        return _vdata.begin();
      }

    //! Get the end iterator.
    std::vector<OBGenericData*>::iterator GetEnd()
      {
        return _vdata.end();
      }

    //! Delete the matching OBGenericData element.
    void DeleteData(OBGenericData *gd)
    {
      std::vector<OBGenericData*>::iterator i;
      for (i = _vdata.begin();i != _vdata.end();++i)
        if (*i == gd)
          {
            delete *i;
            _vdata.erase(i);
            break; // Done, don't do anything more, since iterator is invalid
          }
    }

  }; // OBSetData

  //! \class OBVirtualBond generic.h <openbabel/generic.h>
  //! \brief Used to temporarily store bonds that reference
  //! an atom that has not yet been added to a molecule
 class OBAPI OBVirtualBond : public OBGenericData
  {
  protected:
    int _bgn;
    int _end;
    int _ord;
    int _stereo;
  public:
    OBVirtualBond();
    virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new OBVirtualBond(*this);}
    OBVirtualBond(int,int,int,int stereo=0);
    int GetBgn()    {      return(_bgn);    }
    int GetEnd()    {      return(_end);    }
    int GetOrder()  {      return(_ord);    }
    int GetStereo() {      return(_stereo); }
  };

  //! \class OBRingData generic.h <openbabel/generic.h>
  //! \brief Used to store the SSSR set (filled in by OBMol::GetSSSR())
 class OBAPI OBRingData : public OBGenericData
  {
  protected:
    std::vector<OBRing*> _vr;
  public:
    OBRingData();
    OBRingData(const OBRingData &);
    virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new OBRingData(*this);}
    ~OBRingData();

    OBRingData &operator=(const OBRingData &);

    void SetData(std::vector<OBRing*> &vr)
    {
      _vr = vr;
    }
    void PushBack(OBRing *r)
    {
      _vr.push_back(r);
    }
    std::vector<OBRing*> &GetData()
      {
        return(_vr);
      }

    std::vector<OBRing*>::iterator BeginRings()
      { return(_vr.begin()); }
    std::vector<OBRing*>::iterator EndRings()
      { return(_vr.end()); }
    OBRing *BeginRing(std::vector<OBRing*>::iterator &i);
    OBRing *NextRing(std::vector<OBRing*>::iterator &i);
  };

  //! \class OBUnitCell generic.h <openbabel/generic.h>
  //! \brief Used for storing information about periodic boundary conditions
  //!   with conversion to/from translation vectors and
  //!  (a, b, c, alpha, beta, gamma)
 class OBAPI OBUnitCell: public OBGenericData
  {
  public:
    enum LatticeType { Undefined, 
                       Triclinic, 
                       Monoclinic, 
                       Orthorhombic, 
                       Tetragonal, 
                       Rhombohedral /**< also called trigonal*/, 
                       Hexagonal, 
                       Cubic};


  protected:
    matrix3x3 _mOrtho;// Orthogonal matrix of column vectors 
    matrix3x3 _mOrient;// Orientation matrix
    vector3 _offset;
    std::string _spaceGroupName;
    const SpaceGroup* _spaceGroup;
    LatticeType _lattice;
  public:
    //! public constructor
    OBUnitCell();
    OBUnitCell(const OBUnitCell &);
    virtual OBGenericData* Clone(OBBase* /*parent*/) const
    {return new OBUnitCell(*this);}
    ~OBUnitCell()    {}

    OBUnitCell &operator=(const OBUnitCell &);

    /*!
    **\brief Constructs the cell matrix in lower triangular form from the values supplied
    **\param a The length a
    **\param b The length b
    **\param c The length c
    **\param alpha The angle alpha
    **\param beta The angle beta
    **\param gamma The angle gamma
    */
    void SetData(const double a, const double b, const double c,
                 const double alpha, const double beta, const double gamma);
    /*!
    **\brief Constructs the cell matrix using the supplied row vectors
    **\param v1 The x-vector
    **\param v2 The y-vector
    **\param v3 The z-vector
    **\see OBUnitCell::GetCellVectors
    */
    void SetData(const vector3 v1, const vector3 v2, const vector3 v3);

    /*!
    **\brief Sets the unit cell matrix
    **\param m The unit cell matrix (row vectors)
    **\see OBUnitCell::GetCellMatrix
    */
    void SetData(const matrix3x3 m);

    //! Set the offset to the origin to @p v1
    void SetOffset(const vector3 v1);

    //! Set the space group for this unit cell.
    //! Does not create an OBSymmetryData entry
    void SetSpaceGroup(const SpaceGroup* sg) { _spaceGroup = sg; }

    //! Set the space group symbol for this unit cell.
    //! Does not create an OBSymmetryData entry or attempt to convert
    //!  between different symbol notations
    void SetSpaceGroup(const std::string sg) { _spaceGroup = SpaceGroup::GetSpaceGroup (sg); 
                                               _spaceGroupName = sg; }
    
    //! Set the space group for this unit cell. Each spacegroup-symbol
    //! has a numeric equivalent which is easier to use to convert between
    //! notations.
    //! Does not create an OBSymmetryData entry or attempt to convert
    //!  between different symbol notations
    void SetSpaceGroup(const int sg) { _spaceGroup = SpaceGroup::GetSpaceGroup (sg); }
    
    //! Set the Bravais lattice type for this unit cell
    void SetLatticeType(const LatticeType lt) { _lattice = lt; }

    //! Duplicate symmetry-unique atoms to fill out the unit cell
    //! of the molecule, based on the known space group
    void FillUnitCell(OBMol *);

    //! \return vector a
    double GetA();
    //! \return vector b
    double GetB();
    //! \return vector c
    double GetC();
    //! \return angle alpha
    double GetAlpha();
    //! \return angle beta
    double GetBeta();
    //! \return angle gamma
    double GetGamma();
    //! \return any offset in the origin of the periodic boundaries
    vector3 GetOffset();

    //! \return the text representation of the space group for this unit cell
    const SpaceGroup* GetSpaceGroup() { return(_spaceGroup); }
		
    //! \return the text representation of the space group for this unit cell
    const std::string GetSpaceGroupName() { return(_spaceGroupName); }

    //! \return lattice type (based on the @p spacegroup)
    LatticeType GetLatticeType( int spacegroup );
    
    //! \return lattice type (based on angles and distances)
    LatticeType GetLatticeType();

    //! \return v1, v2, v3 cell vectors
    std::vector<vector3> GetCellVectors();
    //! Access to the cell matrix as row vectors, useful for writing input files.
    //! Equivalent to the transpose of GetOrientationMatrix() * GetOrthoMatrix()
    //! \return The cell matrix with row vectors.
    //! \see OBUnitCell::GetOrthoMatrix
    //! \see OBUnitCell::GetFractionalMatrix
    //! \see OBUnitCell::GetOrientationMatrix
    //! \see OBUnitCell::FractionalToCartesian
    //! \see OBUnitCell::CartesianToFractional
    matrix3x3	GetCellMatrix();
    //! \return The orthogonalization matrix, used for converting from fractional to Cartesian coords.
    //! \see OBUnitCell::GetCellMatrix
    //! \see OBUnitCell::GetFractionalMatrix
    //! \see OBUnitCell::GetOrientationMatrix
    //! \see OBUnitCell::FractionalToCartesian
    //! \see OBUnitCell::CartesianToFractional
    matrix3x3 GetOrthoMatrix();
    //! Used to convert fractional and cartesian coordinates if the
    //! cell is not oriented in standard form (a parallel to x axis,
    //! b in xy plane)
    //! \return The orientation matrix
    //! \see OBUnitCell::GetOrthoMatrix
    //! \see OBUnitCell::GetCellMatrix
    //! \see OBUnitCell::GetFractionalMatrix
    //! \see OBUnitCell::FractionalToCartesian
    //! \see OBUnitCell::CartesianToFractional
    matrix3x3 GetOrientationMatrix();
    //! \return The fractionalization matrix, used for converting from Cartesian to fractional coords.
    //! \see OBUnitCell::GetOrthoMatrix
    //! \see OBUnitCell::GetCellMatrix
    //! \see OBUnitCell::GetOrientationMatrix
    //! \see OBUnitCell::FractionalToCartesian
    //! \see OBUnitCell::CartesianToFractional
    matrix3x3 GetFractionalMatrix();

    //! Convenience function to convert fractional coordinates to 
    //! cartesian coordinates. Returns
    //! 
    //! GetOrientationMatrix() * GetOrthoMatrix() * frac + GetOffset()
    //! \param frac Vector of fractional coordinates
    //! \return Cartesian coordinates
    vector3 FractionalToCartesian(vector3 frac);
    //! Convenience function to convert cartesian coordinates to 
    //! fractional coordinates. Returns
    //! 
    //! GetFractionalMatrix() * GetOrientationMatrix().inverse() * (cart - GetOffset())
    //! \param cart Vector of cartesian coordinates
    //! \return Fractional coordinates
    vector3 CartesianToFractional(vector3 cart);

    //! Wraps cartesian coordinate to fall within the unit cell.
    //! \param cart Vector of cartesian coordinates
    //! \return Cartesian coordinates within cell boundaries.
    vector3 WrapCartesianCoordinate(vector3 cart);
    //! Wraps fractional coordinate to fall within the unit cell.
    //! \param frac Vector of fractional coordinates
    //! \return Fractional coordinates within cell boundaries (between 0 and 1).
    vector3 WrapFractionalCoordinate(vector3 frac);

    //! \return The numeric value of the given spacegroup
    int GetSpaceGroupNumber( std::string name = "" );
    //! \return The cell volume (in Angstroms^3)
    double GetCellVolume();
  };

  //! \class OBConformerData generic.h <openbabel/generic.h>
  //! \brief Used to hold data on conformers or geometry optimization steps
  //!
  //! Supplements the support for multiple coordinate sets in OBMol, e.g.,
  //! OBMol::GetConformer()
 class OBAPI OBConformerData: public OBGenericData
  {
  protected:
    //! Dimensionalities of conformers
    std::vector<unsigned short>              _vDimension;
    //! Relative energies of conformers (preferably in kJ/mol)
    std::vector<double>                      _vEnergies;
    //! Atomic forces for each conformer
    std::vector< std::vector< vector3 > >    _vForces;
    //! Atomic velocities for each conformer (e.g., trajectories)
    std::vector< std::vector< vector3 > >    _vVelocity;
    //! Atomic displacements for each conformer (e.g., RMS distances)
    std::vector< std::vector< vector3 > >    _vDisplace;
    //! Additional data (as strings)
    std::vector<std::string>                 _vData;
    
  public:
    OBConformerData();
    OBConformerData(const OBConformerData &);
    virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new OBConformerData(*this);}
    ~OBConformerData()    {}

    OBConformerData &operator=(const OBConformerData &);

    void SetDimension(std::vector<unsigned short> vd) { _vDimension = vd; }
    void SetEnergies(std::vector<double> ve) { _vEnergies = ve; }
    void SetForces(std::vector< std::vector< vector3 > > vf) {_vForces = vf;}
    void SetVelocities(std::vector< std::vector< vector3 > > vv)
    { _vVelocity = vv; }
    void SetDisplacements(std::vector< std::vector< vector3 > > vd)
    { _vDisplace = vd; }
    void SetData(std::vector<std::string> vdat) { _vData = vdat; }

    std::vector<unsigned short> GetDimension() { return _vDimension; }
    std::vector<double>         GetEnergies()  { return _vEnergies; }
    std::vector< std::vector< vector3 > > GetForces() {return _vForces; }
    std::vector< std::vector< vector3 > > GetVelocities()
      {return _vVelocity;}
    std::vector< std::vector< vector3 > > GetDisplacements()
      {return _vDisplace;}
    std::vector<std::string>    GetData() { return _vData; }

  };

  //! \class OBSymmetryData generic.h <openbabel/generic.h>
  //! \brief Used to hold the point-group and/or space-group symmetry
  //! \todo Add support for translation between symbol notations.
  //!        Add symmetry perception routines.
 class OBAPI OBSymmetryData: public OBGenericData
  {
  protected:
    std::string _spaceGroup;
    std::string _pointGroup;
  public:
    OBSymmetryData();
    OBSymmetryData(const OBSymmetryData &);
    virtual OBGenericData* Clone(OBBase* /*parent*/) const{return new OBSymmetryData(*this);}
    ~OBSymmetryData()    {}

    OBSymmetryData &operator=(const OBSymmetryData &);

    void SetData(std::string pg, std::string sg = "")
    { _pointGroup = pg; _spaceGroup = sg; }
    void SetPointGroup(std::string pg) { _pointGroup = pg; }
    void SetSpaceGroup(std::string sg) { _spaceGroup = sg; }

    std::string GetPointGroup() { return _pointGroup; }
    std::string GetSpaceGroup() { return _spaceGroup; }
  };

  //! \class OBTorsion generic.h <openbabel/generic.h>
  //! \brief Used to hold the torsion data for a single rotatable bond
  //! and all four atoms around it
  class OBAPI OBTorsion
  {
    friend class OBMol;
    friend class OBTorsionData;

  protected:
    std::pair<OBAtom*,OBAtom*> _bc;
    //! double is angle in radians
    std::vector<triple<OBAtom*,OBAtom*,double> > _ads;

    OBTorsion(): _bc(NULL, NULL)      {      }
    //protected for use only by friend classes
    OBTorsion(OBAtom *, OBAtom *, OBAtom *, OBAtom *);

    std::vector<quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> > GetTorsions();

  public:
    OBTorsion(const OBTorsion &);
    ~OBTorsion()      {}

    OBTorsion& operator=(const OBTorsion &);

    void Clear();
    bool Empty()    {      return(_bc.first == 0 && _bc.second == 0);    }

    bool AddTorsion(OBAtom *a,OBAtom *b, OBAtom *c,OBAtom *d);
    bool AddTorsion(quad<OBAtom*,OBAtom*,OBAtom*,OBAtom*> &atoms);

    bool SetAngle(double radians, unsigned int index = 0);
    bool SetData(OBBond * /*bond*/) { return false; }

    bool GetAngle(double &radians, unsigned int index =0);
    //! Gets the bond index of the central bond
    //! \return int bond index
    unsigned int GetBondIdx();
    unsigned int GetSize() const    {      return _ads.size();    }

    //! Gets the two central atoms of ABCD torsion
    //!   \return pair<OBAtom*,OBAtom*>
    std::pair<OBAtom*,OBAtom*>                  GetBC()
      {
        return(_bc);
      }
    //! Gets the vector of distal atoms of ABCD torsion
    //! \return vector of A,D atom pointers and a double
    std::vector<triple<OBAtom*,OBAtom*,double> > GetADs()
    {
      return(_ads) ;
    }

    bool IsProtonRotor();
  };

  //! \class OBTorsionData generic.h <openbabel/generic.h>
  //! \brief Used to hold torsions as generic data for OBMol.
  //! 
  //! Filled by OBMol::FindTorsions()
 class OBAPI OBTorsionData : public OBGenericData
  {
    friend class OBMol;

  protected:
    std::vector<OBTorsion> _torsions;

    OBTorsionData();
    OBTorsionData(const OBTorsionData &);

  public:
    OBTorsionData &operator=(const OBTorsionData &);

    //! \todo Needs to be updated to rebase atom pointers (or use indexes)
    virtual OBGenericData* Clone(OBBase* /*parent*/) const
    {return new OBTorsionData(*this);}

    void Clear();

    //! Gets a vector of the OBTorsion objects
    //! \return the vector of torsions
    std::vector<OBTorsion> GetData() const
      {
        return _torsions;
      }
      
    //! Gets the number of torsion structs
    //! \return integer count of the number of torsions
    unsigned int      GetSize() const
    {
      return _torsions.size();
    }

    void SetData(OBTorsion &torsion);

    bool FillTorsionArray(std::vector<std::vector<unsigned int> > &torsions);
  };

  //! \class OBAngle generic.h <openbabel/generic.h>
  //! \brief Used to hold the 3 atoms in an angle and the angle itself
  class OBAPI OBAngle
  {
    friend class OBMol;
    friend class OBAngleData;

  protected:

    //member data

    OBAtom                *_vertex;
    std::pair<OBAtom*,OBAtom*>  _termini;
    double                  _radians;

    //protected member functions

    OBAngle();	//protect constructor for use only by friend classes
    OBAngle(OBAtom *vertex,OBAtom *a,OBAtom *b);

    triple<OBAtom*,OBAtom*,OBAtom*> GetAtoms();
    void SortByIndex();

  public:

    OBAngle(const OBAngle &);
    ~OBAngle()
      {
        _vertex = NULL;
      }

    OBAngle &operator = (const OBAngle &);
    bool     operator ==(const OBAngle &);

    void  Clear();

    //! Gets the OBAngle angle value
    //! \return angle in radians
    double GetAngle() const
    {
      return(_radians);
    }
    //! Sets the OBAngle to @p radians
    //! \param angle in radians
    void  SetAngle(double angle)
    {
      _radians = angle;
    }
    void  SetAtoms(OBAtom *vertex,OBAtom *a,OBAtom *b);
    void  SetAtoms(triple<OBAtom*,OBAtom*,OBAtom*> &atoms);

  };

  //! \class OBAngleData generic.h <openbabel/generic.h>
  //! \brief Used to hold all angles in a molecule as generic data for OBMol
 class OBAPI OBAngleData : public OBGenericData
  {
    friend class OBMol;

  protected:
    std::vector<OBAngle> _angles;

    OBAngleData();
    OBAngleData(const OBAngleData &);
    //! Gets the angle vector data
    /** \return a vector<OBAngle> **/
    std::vector<OBAngle> GetData() const
      {
        return(_angles);
      }

  public:
    OBAngleData &operator =(const OBAngleData &);
    virtual OBGenericData* Clone(OBBase* /*parent*/) const
    {return new OBAngleData(*this);}

    void Clear();
    unsigned int FillAngleArray(int **angles, unsigned int &size);
    bool FillAngleArray(std::vector<std::vector<unsigned int> > &angles);

    void         SetData(OBAngle &);
    //! Gets the number of angles stored
    //! \return integer count of the number of angles
    unsigned int GetSize() const
    {
      return _angles.size();
    }
  };

  enum atomreftype{
    output,     //!< 
    input,      //!<
    calcvolume  //!<
  }; // sets which atom4ref is accessed by OBChiralData

  //! \class OBChiralData generic.h <openbabel/generic.h>
  //! \brief Used to hold chiral inforamtion about the atom as OBGenericData
 class OBAPI OBChiralData : public OBGenericData
  {
    friend class OBMol;
    friend class OBAtom;

  protected:
    std::vector<unsigned int> _atom4refs; //!< input atom references
    std::vector<unsigned int> _atom4refo; //!< output atom references
    std::vector<unsigned int> _atom4refc; //!< calcvolume references

    //! The parity of the vector (of length 4)
    //! 1234 returns 0, 1243 returns 1
    int parity;

  public:

    OBChiralData();
    OBChiralData(const OBChiralData &src);
    virtual OBGenericData* Clone(OBBase* /*parent*/) const
      { return new OBChiralData(*this); }
    OBChiralData &operator =(const OBChiralData &);
    ~OBChiralData(){}

    void Clear();

    //! \return a vector of all 4 atom references of type @p t
    std::vector<unsigned int> GetAtom4Refs(atomreftype t) const;
    //! \return the atom reference specified by @p a of type @p t
    unsigned int GetAtomRef(int a,atomreftype t);

    bool SetAtom4Refs(std::vector<unsigned int> atom4refs, atomreftype t);
    int AddAtomRef(unsigned int atomref, atomreftype t);

    //! \return the size of the atom references of a given type @p t
    unsigned int GetSize(atomreftype t) const;
  };

  //! \class OBSerialNums generic.h <openbabel/generic.h>
  //! \brief Defines a map between serial numbers (e.g., in a PDB file) and OBAtom objects inside a molecule
 class OBSerialNums : public OBGenericData
  {
  protected:
    std::map<int, OBAtom*> _serialMap; //!< map between serial num

  public:

  OBSerialNums() :
    OBGenericData("obSerialNums", OBGenericDataType::SerialNums)
      {}

  OBSerialNums(const OBSerialNums &cp) : OBGenericData(cp)
    {
      _serialMap = cp._serialMap;
    }
    //! Member variables contain OBAtom pointers, so copying only valid within
    //! same molecule, unless the code is modified, as in OBRotamerList 
    virtual OBGenericData* Clone(OBBase* /*parent*/) const
    {return new OBSerialNums(*this);}
      
    std::map<int,OBAtom*> &GetData()    { return _serialMap;    }
    void SetData(std::map<int,OBAtom*> &sm) { _serialMap = sm;  }
      
  };

  //! \class OBVibrationData generic.h <openbabel/generic.h>
  //! \brief Used to hold the normal modes of a molecule, etc.
 class OBAPI OBVibrationData: public OBGenericData
  {
  protected:
    //! Normal modes in 1/sqrt(a.u.)
    std::vector< std::vector< vector3 > > _vLx;
    
    //! Harmonic frequencies in inverse centimeters
    std::vector<double>  _vFrequencies;
    
    //! Infrared absorption intensities in KM/Mole
    std::vector<double>  _vIntensities;

    //! Raman activities  
    std::vector<double>  _vRamanActivities;

  public:
    OBVibrationData(): OBGenericData("VibrationData", OBGenericDataType::VibrationData){};
    virtual ~OBVibrationData() {}
    virtual OBGenericData* Clone(OBBase*) const
         {return new OBVibrationData(*this);}
    
    OBVibrationData & operator=(const OBVibrationData &);
    
    void SetData(const std::vector< std::vector< vector3 > > & lx,
                 const std::vector<double> & frequencies,
                 const std::vector<double> & intensities);
    void SetData(const std::vector< std::vector< vector3 > > &,
                 const std::vector<double> &,
                 const std::vector<double> &,
                 const std::vector<double> &);
    
    std::vector< std::vector< vector3 > > GetLx() const
      { return this->_vLx; }
    std::vector<double> GetFrequencies() const
      { return this->_vFrequencies; }
    std::vector<double> GetIntensities() const
      { return this->_vIntensities; }
    std::vector<double> GetRamanActivities() const
      { return this->_vRamanActivities; }

    unsigned int GetNumberOfFrequencies() const;
};

  //! \class OBDOSData generic.h <openbabel/generic.h>
  //! \brief Used to hold density of states information
 class OBAPI OBDOSData: public OBGenericData
  {
  protected:
    //! Fermi energy (eV) as shown in _vEnergies
    double _fermi;
    
    //! Energy levels (eV)
    std::vector<double>  _vEnergies;
    
    //! Density of corresponding energy level (number of states / unit cell)
    std::vector<double>  _vDensities;

    //! Integrated DOS vector
    std::vector<double>  _vIntegration;

  public:
    OBDOSData(): OBGenericData("DOSData", OBGenericDataType::DOSData){};
    virtual ~OBDOSData() {}
    virtual OBGenericData* Clone(OBBase*) const
         {return new OBDOSData(*this);}
    
    OBDOSData & operator=(const OBDOSData &);
    
    void SetData(double,
                 const std::vector<double> &,
                 const std::vector<double> &,
                 const std::vector<double> &);
    
    double GetFermiEnergy() const
      { return this->_fermi; }
    std::vector<double> GetEnergies() const
      { return this->_vEnergies; }
    std::vector<double> GetDensities() const
      { return this->_vDensities; }
    std::vector<double> GetIntegration() const
      { return this->_vIntegration; }
};

  //! \class OBOrbitalEnergyData generic.h <openbabel/generic.h>
  //! \brief Used to hold information about orbital energies
  class OBAPI OBOrbitalEnergyData: public OBGenericData
  {
  protected:
    std::vector<double> _alphaEigenvalues;
    std::vector<double> _betaEigenvalues;
    std::vector<std::string> _alphaMullikenSymbols;
    std::vector<std::string> _betaMullikenSymbols;
  };

  //! \class OBElectronicTransitionData generic.h <openbabel/generic.h>
  //! \brief Used to hold information about electronic transitions
 class OBAPI OBElectronicTransitionData: public OBGenericData
  {
  protected:
    //! Wavelengths (nm)
    std::vector<double>  _vWavelengths;
    
    //! Oscillator strengths
    std::vector<double>  _vForces;

    //! Electric dipole strengths
    std::vector<double>  _vEDipole;

    //! Rotatory strengths (velocity)
    std::vector<double>  _vRotatoryStrengthsVelocity;

    //! Rotatory strengths (length)
    std::vector<double>  _vRotatoryStrengthsLength;

  public:
    OBElectronicTransitionData(): OBGenericData("ElectronicTransitionData", OBGenericDataType::ElectronicTransitionData) {}
    virtual ~OBElectronicTransitionData() {}
    virtual OBGenericData* Clone(OBBase*) const
         {return new OBElectronicTransitionData(*this);}
    
    OBElectronicTransitionData & operator=(const OBElectronicTransitionData &);
    
    void SetData(const std::vector<double> & wavelengths,
                 const std::vector<double> & forces);

    void SetEDipole(const std::vector<double> &);
    void SetRotatoryStrengthsVelocity(const std::vector<double> &);
    void SetRotatoryStrengthsLength(const std::vector<double> &);

    std::vector<double> GetWavelengths() const
      { return this->_vWavelengths; }
    std::vector<double> GetForces() const
      { return this->_vForces; }
    std::vector<double> GetEDipole() const
      { return this->_vEDipole; }
    std::vector<double> GetRotatoryStrengthsVelocity() const
      { return this->_vRotatoryStrengthsVelocity; }
    std::vector<double> GetRotatoryStrengthsLength() const
      { return this->_vRotatoryStrengthsLength; }
};

  //! \class OBRotationData generic.h <openbabel/generic.h>
  //! \brief Used to hold the rotational constants and symmetry numbers
 class OBAPI OBRotationData: public OBGenericData
 {
 public:
   enum RType{UNKNOWN, ASYMMETRIC, SYMMETRIC, LINEAR};
   OBRotationData(): OBGenericData("RotationData", OBGenericDataType::RotationData){}
   virtual ~OBRotationData(){};
   virtual OBGenericData* Clone(OBBase*) const
         {return new OBRotationData(*this);}
   void SetData(RType RotorType, std::vector<double> RotationalConstants, int SymmetryNumber)
   {
     RotConsts = RotationalConstants;
     type = RotorType;
     SymNum = SymmetryNumber;
   }

   /// \return Rotational constants in GHz
   std::vector<double> GetRotConsts()const{ return RotConsts; }

   int GetSymmetryNumber()const{ return SymNum; }
   RType GetRotorType()const   { return type; }

 protected:
   std::vector<double> RotConsts;//!< Rotational constants in GHz
   int                 SymNum;   //!< Rotational Symmetry Number
   RType               type;     //!< linear, symmetric or asymmetric top
 };
 
  //! \class OBVectorData generic.h <openbabel/generic.h>
  //! \brief Used to hold a 3D vector item (e.g., a dipole moment)
  //! \since version 2.2
 class OBAPI OBVectorData: public OBGenericData
 {
 public:
   OBVectorData(): OBGenericData("VectorData", OBGenericDataType::VectorData){}
   virtual ~OBVectorData(){};
   virtual OBGenericData* Clone(OBBase*) const
         {return new OBVectorData(*this);}
   void SetData(double x, double y, double z)
     { _vec = vector3(x, y, z); }
   void SetData(vector3 data)
     { _vec = data; }
   vector3 GetData() const
     { return _vec; }
   
 protected:
   vector3            _vec; //!< 3D vector to be stored
 };
 
   //! \class OBMatrixData generic.h <openbabel/generic.h>
   //! \brief Used to hold a 3x3 matrix item (e.g., a quadrupole moment)
   //! \since version 2.2
  class OBAPI OBMatrixData: public OBGenericData
  {
  public:
    OBMatrixData(): OBGenericData("MatrixData", OBGenericDataType::MatrixData){}
    virtual ~OBMatrixData(){};
    virtual OBGenericData* Clone(OBBase*) const
          {return new OBMatrixData(*this);}
    void SetData(matrix3x3 data)
      { _matrix = data; }
    matrix3x3 GetData() const
      { return _matrix; }

  protected:
    matrix3x3            _matrix; //!< 3x3 matrix to be stored
  };

 //! A standard iterator over vectors of OBGenericData (e.g., inherited from OBBase)
  typedef std::vector<OBGenericData*>::iterator OBDataIterator;

} //end namespace OpenBabel

#endif // OB_GENERIC_H

//! \file generic.h
//! \brief Handle generic data classes. Custom data for atoms, bonds, etc.