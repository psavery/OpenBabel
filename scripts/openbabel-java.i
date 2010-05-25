%module openbabel_java

%{
// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif


#include <openbabel/obutil.h>
#include <openbabel/rand.h>
#include <openbabel/math/vector3.h>
#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
#include <openbabel/math/spacegroup.h>

#include <openbabel/generic.h>
#include <openbabel/griddata.h>

#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/residue.h>
#include <openbabel/internalcoord.h>

#include <openbabel/ring.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/plugin.h>
#include <openbabel/fingerprint.h>
#include <openbabel/descriptor.h>
#include <openbabel/format.h>

#include <openbabel/forcefield.h>
#include <openbabel/builder.h>
#include <openbabel/op.h>

#include <openbabel/bitvec.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <openbabel/alias.h>
#include <openbabel/atomclass.h>

#include <openbabel/kinetics.h>
#include <openbabel/rotor.h>
#include <openbabel/rotamer.h>

%}

%include "std_map.i"
%include "std_vector.i"
%include "std_string.i"

namespace std {

%define VVTEMPLATE_WRAP(name, T) 
%feature("ignore") vector< vector<T> >::append;
%feature("ignore") vector< vector<T> >::assign;
%feature("ignore") vector< vector<T> >::back;
%feature("ignore") vector< vector<T> >::begin;
%feature("ignore") vector< vector<T> >::capacity;
%feature("ignore") vector< vector<T> >::clear;
%feature("ignore") vector< vector<T> >::empty;
%feature("ignore") vector< vector<T> >::end;
%feature("ignore") vector< vector<T> >::erase;
%feature("ignore") vector< vector<T> >::front;
%feature("ignore") vector< vector<T> >::get_allocator;
%feature("ignore") vector< vector<T> >::insert;
%feature("ignore") vector< vector<T> >::pop;
%feature("ignore") vector< vector<T> >::pop_back;
%feature("ignore") vector< vector<T> >::push_back;
%feature("ignore") vector< vector<T> >::rbegin;
%feature("ignore") vector< vector<T> >::rend;
%feature("ignore") vector< vector<T> >::reserve;
%feature("ignore") vector< vector<T> >::resize;
%feature("ignore") vector< vector<T> >::size;
%feature("ignore") vector< vector<T> >::swap;
%template(vectorv ## name) vector< vector<T> >;
%enddef

%define VECTORTEMPLATE_WRAP(vectorname, T) 
%feature("ignore") vector<T>::append;
%feature("ignore") vector<T>::assign;
%feature("ignore") vector<T>::back;
%feature("ignore") vector<T>::begin;
%feature("ignore") vector<T>::capacity;
%feature("ignore") vector<T>::clear;
%feature("ignore") vector<T>::empty;
%feature("ignore") vector<T>::end;
%feature("ignore") vector<T>::erase;
%feature("ignore") vector<T>::front;
%feature("ignore") vector<T>::get_allocator;
%feature("ignore") vector<T>::insert;
%feature("ignore") vector<T>::pop;
%feature("ignore") vector<T>::pop_back;
%feature("ignore") vector<T>::push_back;
%feature("ignore") vector<T>::rbegin;
%feature("ignore") vector<T>::rend;
%feature("ignore") vector<T>::reserve;
%feature("ignore") vector<T>::resize;
%feature("ignore") vector<T>::size;
%feature("ignore") vector<T>::swap;
%template(vector ## vectorname) vector<T>;
%enddef

VECTORTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(UnsignedInt, unsigned int)
VVTEMPLATE_WRAP(Int, int)
VECTORTEMPLATE_WRAP(Double, double)
VECTORTEMPLATE_WRAP(String, std::string)
VECTORTEMPLATE_WRAP(Vector3, OpenBabel::vector3)
VECTORTEMPLATE_WRAP(OBMol, OpenBabel::OBMol)
VECTORTEMPLATE_WRAP(OBBond, OpenBabel::OBBond)
VECTORTEMPLATE_WRAP(OBResidue, OpenBabel::OBResidue)
VECTORTEMPLATE_WRAP(OBRing, OpenBabel::OBRing)
VECTORTEMPLATE_WRAP(pOBRing, OpenBabel::OBRing*)
VECTORTEMPLATE_WRAP(pOBGenericData, OpenBabel::OBGenericData*)
}

%define CAST_GENERICDATA_TO(subclass)
%inline %{
OpenBabel::OB ## subclass *to ## subclass(OpenBabel::OBGenericData *data) {
    return (OpenBabel::OB ## subclass *) data;
}
%}
%enddef
CAST_GENERICDATA_TO(AngleData)
CAST_GENERICDATA_TO(AtomClassData)
CAST_GENERICDATA_TO(ChiralData)
CAST_GENERICDATA_TO(CommentData)
CAST_GENERICDATA_TO(ConformerData)
CAST_GENERICDATA_TO(ExternalBondData)
CAST_GENERICDATA_TO(GridData)
CAST_GENERICDATA_TO(MatrixData)
CAST_GENERICDATA_TO(NasaThermoData)
CAST_GENERICDATA_TO(PairData)
// CAST_GENERICDATA_TO(PairTemplate)
CAST_GENERICDATA_TO(RateData)
CAST_GENERICDATA_TO(RotamerList)
CAST_GENERICDATA_TO(RotationData)
CAST_GENERICDATA_TO(SerialNums)
CAST_GENERICDATA_TO(SetData)
CAST_GENERICDATA_TO(SymmetryData)
CAST_GENERICDATA_TO(TorsionData)
CAST_GENERICDATA_TO(UnitCell)
CAST_GENERICDATA_TO(VectorData)
CAST_GENERICDATA_TO(VibrationData)
CAST_GENERICDATA_TO(VirtualBond)

%ignore *::operator=;
%ignore *::operator*;
%ignore *::operator*=;
%ignore *::operator+;
%ignore *::operator+=;
%ignore *::operator-;
%ignore *::operator-=;
%ignore *::operator++;
%ignore *::operator--;
%ignore *::operator/;
%ignore *::operator/=;
%ignore *::operator bool;
%ignore *::operator[];

%import <openbabel/babelconfig.h>

%include <openbabel/data.h>
%include <openbabel/rand.h>
%include <openbabel/obutil.h>
%include <openbabel/math/vector3.h>
%warnfilter(503) OpenBabel::matrix3x3; // Not wrapping any of the overloaded operators
%include <openbabel/math/matrix3x3.h>
%include <openbabel/math/transform3d.h>
%include <openbabel/math/spacegroup.h>

%include <openbabel/base.h>
%include <openbabel/generic.h>
%include <openbabel/griddata.h>

%include <openbabel/chains.h>
%include <openbabel/typer.h>

%include <openbabel/plugin.h>

%include <openbabel/oberror.h>
%include <openbabel/format.h>
%include <openbabel/obconversion.h>
%include <openbabel/residue.h>
%include <openbabel/internalcoord.h>

%typemap(javacode) OpenBabel::OBAtom
%{
    private int currentDepth = 0;
    public void SetCurrentDepth(int d) {currentDepth = d;}
    public int GetCurrentDepth() {return currentDepth;}
%}
%include <openbabel/atom.h>
%include <openbabel/bond.h>
%include <openbabel/mol.h>
%include <openbabel/ring.h>
%include <openbabel/parsmart.h>
%include <openbabel/alias.h>
%include <openbabel/atomclass.h>
%ignore OpenBabel::FptIndex;
%include <openbabel/fingerprint.h>
%ignore OpenBabel::OBDescriptor::LessThan;
%include <openbabel/descriptor.h>

# Ignore shadowed methods
%ignore OpenBabel::OBForceField::VectorSubtract(const double *const, const double *const, double *);
%ignore OpenBabel::OBForceField::VectorMultiply(const double *const, const double, double *);
%include <openbabel/forcefield.h>

%include <openbabel/builder.h>
%include <openbabel/op.h>

%warnfilter(503) OpenBabel::OBBitVec; // Not wrapping any of the overloaded operators
%include <openbabel/bitvec.h>

# The following %ignores avoid warning messages due to shadowed classes.
# This does not imply a loss of functionality as (in this case)
# the shadowed class is identical (from the point of view of SWIG) to
# the shadowing class.
# This is because C++ references (&) are transformed by SWIG back into
# pointers, so that OBAtomIter(OBMol &) would be treated the same as
# OBAtomIter(OBMol *).

%ignore OBAtomAtomIter(OBAtom &);
%ignore OBAtomBondIter(OBAtom &);
%ignore OBMolAngleIter(OBMol &);
%ignore OBMolAtomIter(OBMol &);
%ignore OBMolAtomBFSIter(OBMol &, int);
%ignore OBMolAtomBFSIter(OBMol &);
%ignore OBMolAtomDFSIter(OBMol &, int);
%ignore OBMolAtomDFSIter(OBMol &);
%ignore OBMolBondIter(OBMol &);
%ignore OBMolPairIter(OBMol &);
%ignore OBMolRingIter(OBMol &);
%ignore OBMolTorsionIter(OBMol &);
%ignore OBResidueIter(OBMol &);
%ignore OBResidueAtomIter(OBResidue &);

%ignore SpaceGroup::RegisterSpaceGroup;

%define WRAPITERATOR(NAME, RETURNS, OPTIONAL)
 %rename(_next) OpenBabel::NAME::next;
 %rename(inc) OpenBabel::NAME::operator++;
 %rename(HasMore) OpenBabel::NAME::operator bool;
 %typemap(javainterfaces) OpenBabel::NAME "Iterable<RETURNS>, Iterator<RETURNS>";
 %typemap(javaimports) OpenBabel::NAME "  import org.openbabel.RETURNS;
  import java.util.Iterator;"
 %typemap(javacode) OpenBabel::NAME
 %{
	public Iterator<RETURNS> iterator() {
		return this;
	}

	public boolean hasNext() {
		return HasMore();
	}

	public RETURNS next() {
		RETURNS atom = __ref__();
		OPTIONAL
		inc();
		return atom;
	}
	
	public void remove() {
	}	
 %}
%enddef

WRAPITERATOR(OBMolAtomIter, OBAtom, );
WRAPITERATOR(OBMolAtomDFSIter, OBAtom, );
WRAPITERATOR(OBMolAtomBFSIter, OBAtom, atom.SetCurrentDepth(CurrentDepth()););
WRAPITERATOR(OBMolBondIter, OBBond, );
WRAPITERATOR(OBMolAngleIter, vectorUnsignedInt, );
WRAPITERATOR(OBAtomAtomIter, OBAtom, )
WRAPITERATOR(OBAtomBondIter, OBBond, );
WRAPITERATOR(OBMolRingIter, OBRing, );
WRAPITERATOR(OBMolTorsionIter, vectorUnsignedInt, );
WRAPITERATOR(OBResidueIter, OBResidue, );
WRAPITERATOR(OBResidueAtomIter, OBAtom, );
WRAPITERATOR(OBMolPairIter, vectorUnsignedInt, );

%include <openbabel/obiter.h>
