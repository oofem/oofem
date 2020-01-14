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
#ifndef BASICLSRK3D_H
#define BASICLSRK3D_H


#include "../sm/Elements/structuralelement.h"
#include "dofmanager.h"

///@name Input fields for BasicLSR3d
//@{
#define _IFT_BasicLSRK3d_Name "basicLSRK3d"
#define _IFT_BasicLSRK3d_dofstocondense "dofstocondense"
#define _IFT_BasicLSRK3d_refnode "refnode"
#define _IFT_BasicLSRK3d_refangle "refangle"
#define _IFT_BasicLSRK3d_zaxis "zaxis"
//@}

namespace oofem {

class FEI2dQuadLin;

/**
 * This class implements a 2-dimensional basic linear strain rectangle element
 * that can be arbitrary oriented in space, in contrast to basicLSR element
 * that is defined in xy-plane
 * To allow this, third degree of freedom, w-displacement is introduced
 * It is only used for the space orientation of element, because all computations regarding the element
 * are done in local coordinate system, where the element posseses only 2 degrees of freedom
 * This elements also possesses extra stiffness added to its stiffness matrix
 * It is computed by dividing maximum value of stiffness on the diagonal of the stiffness matrix with a factor
 * This allows solving problems of structures with elements in different planes
 * Element is derived using direct integration, numerical integration is not used,
 * stiffness, strain-displacement matrix and displacement interpolation matrix are computed directly
 * from the derived expressions
 * Gauss points used in the element are not used for numerical integration
 * They are used as a reference for computation of stresses
 * crossSectioncharacteristics method returns characteristics of cross section
 * used for calculation of stiffness matrix
 * computeLength method is used for calculation of edge lengths of an element
 * computeLCS method is used for calculation of local coordinate system of an element
 */

class BasicLSRK3d : public StructuralElement
{
protected:
    /// Geometry interpolator only.
    static FEI2dQuadLin interp;

    double a,b,E,ni,t,K;
    int referenceNode;
    FloatMatrix lcsMatrix;
    std::vector< FloatArray > lnodes;



public:
    BasicLSRK3d(int n, Domain *d);
    virtual ~BasicLSRK3d();

    virtual FEInterpolation *giveInterpolation() const;

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);




    virtual int testElementExtension(ElementExtension ext) {
        return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
    }
    //int hasLayeredSupport () {return 1;}

    virtual int computeNumberOfDofs() { return 12; }
    virtual int computeNumberOfGlobalDofs() { return 12; }


    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    // definition & identification
    virtual const char *giveClassName() const { return "BasicLSR3d"; }
    virtual const char *giveInputRecordName() const { return _IFT_BasicLSRK3d_Name; }
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
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &);

    void crossSectioncharacteristics();
    virtual double computeLength();
    virtual void computeClampedStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode rMode, TimeStep *tStep);
    virtual void computeLocalStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode rMode, TimeStep *tStep);
    virtual void computeGaussPoints();

    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 4; }

};
} // end namespace oofem


#endif // BASICLSRK3D_H

