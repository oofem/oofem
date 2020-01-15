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

#ifndef ACM_PLATE_H
#define ACM_PLATE_H



#include "../sm/Elements/structuralelement.h"

#define _IFT_ACMPlate_Name "acmplate"

namespace oofem {
class FEI2dQuadLin;

/**
 * This class implements an quad Discrete Kirchhoff Theory (DKT) element.
 * This element is a plate element suitable for thin plates
 * Element is taken from the book Metoda konačnih elemenata
 * The element has 12 DOFs (w-displacement and rotations along coordinate axes in each node)
 * Element uses direct integration instead of numerical integration
 * Stiffness matrix is derived and the final expression is used for computation
 * Currently element can only compute displacements
 * crossSectioncharacteristics method gives cross section carachtersistics used for computation of stiffness matrix
 * computeLength method gives edge lenghts of an element
 * computeLCS method gives local coordinate system of an element
 */
class ACMPlate : public StructuralElement
{
protected:
    /// Element geometry approximation
    static FEI2dQuadLin interp_lin;
    double a,b,E,ni,t;
    FloatMatrix lcsMatrix;
    std::vector< FloatArray > lnodes;

public:
    ACMPlate(int n, Domain * d);
    virtual ~ACMPlate() { }

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    virtual MaterialMode giveMaterialMode()  { return _2dPlate; }
    virtual int testElementExtension(ElementExtension ext) { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }
	
	// This was put here only to implement the abstract method so that the project could build. The final implementation should be taken care of later. 
	virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) {};
	virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) {};

protected:
    virtual void computeGaussPoints();
	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);
    virtual void computeStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep);

    virtual void giveNodeCoordinates(double &x1, double &x2, double &x3, double &x4,
                                     double &y1, double &y2, double &y3, double &y4,
                                     double &z1, double &z2, double &z3, double &z4);

    void crossSectioncharacteristics();
	virtual double computeLength();

public:
    // definition & identification
    virtual const char *giveClassName() const { return "ACMPlate"; }
    virtual const char *giveInputRecordName() const { return _IFT_ACMPlate_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs() { return 12; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual void computeLCS();
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    int computeLoadGToLRotationMtrx(FloatMatrix &answer);


    virtual void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp);


    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};
} // end namespace oofem


#endif // ACM_PLATE_H

