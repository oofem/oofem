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
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef Tr2Shell7_h
#define Tr2Shell7_h

#include "../sm/Elements/Shells/shell7base.h"
#include "../sm/CrossSections/layeredcrosssection.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_Tr2Shell7_Name "tr2shell7"

namespace oofem {
class FEI3dTrQuad;
class BoundaryLoad;

/**
 * This class represent a 7 parameter shell element.
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */

class Tr2Shell7 : public Shell7Base
{
protected:
    static FEI3dTrQuad interpolation; 
    static IntArray orderingDofTypes;
    static IntArray orderingNodes;
    static IntArray orderingEdgeNodes;

    virtual const IntArray &giveOrderingDofTypes() const;
    virtual const IntArray &giveOrderingNodes() const;
    virtual const IntArray &giveOrderingEdgeNodes() const;

    //specific
    void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;

    virtual double computeVolumeAroundLayer(GaussPoint *mastergp, int layer);
    virtual double computeAreaAround(GaussPoint *gp, double xi);

    virtual void computeGaussPoints();
    //only used for debuging
    void compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer);

    virtual FEInterpolation *giveInterpolation() const;


public:
    Tr2Shell7(int n, Domain * d);        // constructor
    virtual ~Tr2Shell7() { }            // destructor -> declaring as virtual will make each subclass call their respective destr.
    // definition & identification
    virtual int giveNumberOfDofs() { return 42; }
    virtual int giveNumberOfEdgeDofs() { return 21; }
    virtual int giveNumberOfEdgeDofManagers() { return 3;  }
    virtual const char *giveInputRecordName() const { return _IFT_Tr2Shell7_Name; }
    virtual const char *giveClassName() const { return "Tr2Shell7"; }
    //virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_Composite; }
    virtual integrationDomain giveIntegrationDomain() const { return _Triangle; }     // write new wedge-like type 'layeredWedge'
};
} // end namespace oofem
#endif
