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

#ifndef QCSHELL_H
#define QCSHELL_H


#include "../sm/Elements/structuralelement.h"
#include "dofmanager.h"

///@name Input fields for QCSHELL
//@{
#define _IFT_qcshell_Name "qcshell"
#define _IFT_qcshell_dofstocondense "dofstocondense"
#define _IFT_qcshell_refnode "refnode"
#define _IFT_qcshell_refangle "refangle"
#define _IFT_qcshell_zaxis "zaxis"
//@}

namespace oofem {

class FEI2dQuadLin;

/**
 * This class implements an quadrilateral four-node shell element, using Kirchoff's plate theory.
 * Each node has 6 degrees of freedom
 * Element is a combination of CSSRK3d element and ACMPlate element.
 * Its stiffness matrix is a combination of stiffness matrices of these two elements
 * Stiffness related to DOFs 1 and 2 (u,v-displacement) are taken from the membrane element stiffness matrix
 * Stiffness related to DOFS 3,4,5 (w-displacement, u,v-rotations) are taken from the plate element stiffness matrix
 * Stiffness related to DOF 6 (w-rotation) is given as extra stiffness
 * It is computed by dividing maximum stiffness on the stiffness matrix diagonal with a factor
 * Currently element can only compute displacements
 * computeLength method gives edge lengths of an element
 * computeLCS method gives local coordinate system of an element
 */

class qcshell : public StructuralElement
{
protected:
    /// Geometry interpolator only.
    static FEI2dQuadLin interp;

    double a,b,E,ni,t,K;
    int referenceNode;
    FloatMatrix lcsMatrix;
    std::vector< FloatArray > lnodes;



public:
    qcshell(int n, Domain *d);
    virtual ~qcshell();

    virtual FEInterpolation *giveInterpolation() const;

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);




    virtual int testElementExtension(ElementExtension ext) {
        return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
    }
    //int hasLayeredSupport () {return 1;}

    virtual int computeNumberOfDofs() { return 24; }
    virtual int computeNumberOfGlobalDofs() { return 24; }


    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    // definition & identification
    virtual const char *giveClassName() const { return "qcshell"; }
    virtual const char *giveInputRecordName() const { return _IFT_qcshell_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    ///@todo Introduce interpolator and remove these two:
    virtual integrationDomain giveIntegrationDomain() const { return _Square; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_quad_1; }
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    virtual void computeLCS();
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    int computeLoadGToLRotationMtrx(FloatMatrix &answer);


protected:
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);

    void crossSectioncharacteristics();
    virtual double computeLength();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    virtual void computeGaussPoints();

    virtual MaterialMode giveMaterialMode() { return _3dShell; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }

};
} // end namespace oofem


#endif // QCSHELL_H

