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
#ifndef BASICLSR_H
#define BASICLSR_H


#include "../sm/Elements/structuralelement.h"
#include "dofmanager.h"

///@name Input fields for BasicLSR
//@{
#define _IFT_BasicLSR_Name "basicLSR"
#define _IFT_BasicLSR_dofstocondense "dofstocondense"
#define _IFT_BasicLSR_refnode "refnode"
#define _IFT_BasicLSR_refangle "refangle"
#define _IFT_BasicLSR_zaxis "zaxis"
//@}

namespace oofem {

class FEI2dQuadLin;

/**
 * This class implements a 2-dimensional basic linear strain rectangle element
 * Element is taken from the book Ship Structural Design and Analysis
 * Element is applicable for linear-static analysis
 * Element can only be used in x-y plane
 * Element is derived using direct integration, numerical integration is not used,
 * stiffness, strain-displacement matrix and displacement interpolation matrix are computed directly
 * from the derived expressions
 * Gauss points used in the element are not used for numerical integration
 * They are used as a reference for computation of stress
 * crossSectioncharacteristics method returns characteristics of cross section
 * used for calculation of stiffness matrix
 * computeLength method is used for calculation of edge lengths of an element
 * computeLCS method is used for calculation of local coordinate system of an element
 */

class BasicLSR : public StructuralElement
{
protected:
    /// Geometry interpolator only.
    static FEI2dQuadLin interp;

    double a,b,E,ni,t;
    int referenceNode;
    FloatMatrix lcsMatrix;


public:
    BasicLSR(int n, Domain *d);
    virtual ~BasicLSR();

    virtual FEInterpolation *giveInterpolation() const;

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);




    virtual int testElementExtension(ElementExtension ext) {
        return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
    }
    //int hasLayeredSupport () {return 1;}

    virtual int computeNumberOfDofs() { return 8; }

    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    // definition & identification
    virtual const char *giveClassName() const { return "BasicLSR"; }
    virtual const char *giveInputRecordName() const { return _IFT_BasicLSR_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    ///@todo Introduce interpolator and remove these two:
    virtual integrationDomain giveIntegrationDomain() const { return _Square; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_quad_1; }
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    virtual void computeLCS();
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);


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

#endif // BASICLSR_H

