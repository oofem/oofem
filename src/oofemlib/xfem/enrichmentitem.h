/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef enrichmentitem_h
#define enrichmentitem_h

#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"
#include "layeredcrosssection.h"
#include "dofiditem.h"

///@name Input fields for XFEM
//@{
#define _IFT_CrackTip_Name "cracktip"
#define _IFT_CrackInterior_Name "crackinterior"

#define _IFT_Inclusion_Name "inclusion"
#define _IFT_Inclusion_material "material"

#define _IFT_EnrichmentItem_domains "enrichmentdomains"
#define _IFT_EnrichmentItem_function "enrichmentfunction"

#define _IFT_Delamination_Name "delamination"
#define _IFT_Delamination_xiCoords "delaminationxicoords"
//#define _IFT_MultipleDelamination_Name "multipledelamination"
//@}

#define _IFT_Crack_Name "crack"

namespace oofem {
template< class T > class AList;
class BasicGeometry;
class EnrichmentFunction;
class EnrichmentDomain;

/**
 * Abstract class representing entity, which is included in the FE model using one (or more)
 * global functions. Such entity may represent crack, material interface, etc.
 * As the geometry of such entity may be represented in a number of ways, the hierarchy of classes
 * derived from base Geometry class is used to achieve flexibility of geometry representation.
 *
 * Each EnrichmentItem keeps its DOF labels (assigned/allocated by XFemManager, its geometry representation, and
 * keeps the list of its EnrichmentFunctions.
 * @author chamrova
 * @author Jim Brouzoulis
 */
class EnrichmentItem : public FEMComponent
{
public:
    /// Constructor.
    EnrichmentItem(int n, XfemManager *xm, Domain *aDomain);
    virtual ~EnrichmentItem();
    virtual IRResultType initializeFrom(InputRecord *ir);
    int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const = 0;
    IntArray *giveEnrichesDofsWithIdArray() { return this->enrichesDofsWithIdArray; }
    int giveNumberOfEnrDofs();

    // Enrichment domains
    BasicGeometry *giveGeometry(int i);
    BasicGeometry *giveGeometry();
    EnrichmentDomain *giveEnrichmentDomain(int i) { return this->enrichmentDomainList->at(i); }
    int giveNumberOfEnrichmentDomains() { return this->numberOfEnrichmentDomains; }

    // Enrichment functions
    EnrichmentFunction *giveEnrichmentFunction(int n);
    int giveNumberOfEnrichmentfunctions() { return this->numberOfEnrichmentFunctions; }

    // Spatial query
    bool isDofManEnriched(DofManager *dMan);
    bool isDofManEnrichedByEnrichmentDomain(DofManager *dMan, int edNumber);
    bool isElementEnriched(const Element *element); 
    bool isElementEnrichedByEnrichmentDomain(const Element *element, int edNumber); 

    // Should update receiver geometry to the state reached at given time step.
    virtual void updateGeometry(TimeStep *tStep) {};

    int giveStartOfDofIdPool() { return this->startOfDofIdPool; };
    void computeDofManDofIdArray(IntArray &DofIdArray, DofManager *dMan, int enrichmentDomainNumber); // list of id's a particular dof manager supports
    void giveEIDofIdArray(IntArray &answer, int enrichmentDomainNumber); // list of id's for the enrichment dofs


protected:
    /// Link to associated Xfem manager.
    XfemManager *xMan;
    int startOfDofIdPool; // points to the first available dofId number associated with the ei 

    /// Geometry associated with EnrichmentItem.
    IntArray enrichmentDomainNumbers;
    IntArray *enrichesDofsWithIdArray;

    /// EnrichmentFunction associated with the EnrichmentItem. - should generally be a list of functions
    int enrichmentFunction;

    /// Geometry object
    AList< EnrichmentDomain > *enrichmentDomainList;
    int numberOfEnrichmentDomains;

    /// Enrichment function list.
    AList< EnrichmentFunction > *enrichmentFunctionList;
    int numberOfEnrichmentFunctions;

};

/** Sub classes to EnrichmentItem. */
class CrackTip : public EnrichmentItem // only for 2D. Only the tip element belong to this
{
public:
    CrackTip(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain) { }
    virtual const char *giveClassName() const { return "CrackTip"; }
    virtual const char *giveInputRecordName() const { return _IFT_CrackTip_Name; }
};

/** Concrete representation of EnrichmentItem. */
class CrackInterior : public EnrichmentItem // rest of the crack el. that does not contain any tip 
{
public:
    CrackInterior(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain) { }
    virtual const char *giveClassName() const { return "CrackInterior"; }
    virtual const char *giveInputRecordName() const { return _IFT_CrackInterior_Name; }
};

/** Concrete representation of EnrichmentItem. */
class Inclusion : public EnrichmentItem
{
protected:
    Material *mat;
public:
    Inclusion(int n, XfemManager *xm, Domain *aDomain);
    virtual const char *giveClassName() const { return "Inclusion"; }
    virtual const char *giveInputRecordName() const { return _IFT_Inclusion_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Material *giveMaterial() { return mat; }
};


/** Concrete representation of Delamination. */
class Delamination : public EnrichmentItem 
{
public:
    Delamination(int n, XfemManager *xm, Domain *aDomain);
    virtual const char *giveClassName() const { return "Delamination"; }
    virtual const char *giveInputRecordName() const { return _IFT_Delamination_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    FloatArray enrichmentDomainXiCoords; 
    std::list<std::pair<int, double> > delaminationXiCoordList;
    double giveDelaminationZCoord(int n, Element *element); 

    int giveDelaminationGroupAt(double z);
    FloatArray delaminationGroupMidZ(int dGroup);
    double giveDelaminationGroupMidZ(int dGroup, Element *e);
    
    FloatArray delaimnationGroupThickness;
    double giveDelaminationGroupThickness(int dGroup, Element *e);

    void giveDelaminationGroupZLimits(int &dGroup, double &zTop, double &zBottom, Element *e);
    double heaviside(double xi, double xi0);

};

/** Concrete representation of Crack. */
class Crack : public EnrichmentItem
{
//protected:

public:
    Crack(int n, XfemManager *xm, Domain *aDomain);
    virtual const char *giveClassName() const { return "Crack"; }
    virtual const char *giveInputRecordName() const { return _IFT_Crack_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
};

} // end namespace oofem

#endif  // enrichmentitem_h
