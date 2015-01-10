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

#ifndef PrescribedGenStrainShell7_h
#define PrescribedGenStrainShell7_h

#include "boundarycondition.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for PrescribedTensor
//@{
#define _IFT_PrescribedGenStrainShell7_Name "PrescribedGenStrainShell7"
#define _IFT_PrescribedGenStrainShell7_centercoords "ccoord"
#define _IFT_PrescribedGenStrainShell7_generalizedstrain "generalizedstrain"
#define _IFT_PrescribedGenStrainShell7_initialgeneralizedstrain "initialgeneralizedstrain"
//@}

namespace oofem {
/**
 * Prescribes @f$ u_i = d_{ij}(x_j-\bar{x}_j) @f$ or @f$ s = d_{1j}(x_j - \bar{x}_j) @f$
 * where @f$ v_i @f$ are primary unknowns for the subscale.
 * This is typical boundary condition in multiscale analysis where @f$ d = \partial_x s@f$
 * would a macroscopic gradient at the integration point, i.e. this is a boundary condition for prolongation.
 * It is also convenient to use when one wants to test a arbitrary specimen for shear.
 * @author Jim Brouzoulis
 */
class OOFEM_EXPORT PrescribedGenStrainShell7 : public BoundaryCondition
{
protected:
    /// Prescribed gradient @f$ d_{ij} @f$
    FloatMatrix gradient;

    /**
     * Initial generalized strain
     * eps = [Xbar_{,alpha}, M_{,alpha}, M  ] = [ (3  3)  (3  3)  3] - 15 components
     */
    FloatArray initialGenEps;
    /** 
     * Current generalized strain
     * eps = [xbar_{,alpha}, m_{,alpha}, m, gamma_{,alpha}, gamma ] = [ (3  3)  (3  3)  3  (1 1) 1] - 18 components
     */
    FloatArray genEps;

    /// Center coordinate @f$ \bar{x}_i @f$
    FloatArray centerCoord;

public:
    /**
     * Creates boundary condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    PrescribedGenStrainShell7(int n, Domain *d) : BoundaryCondition(n, d) { }

    /// Destructor
    virtual ~PrescribedGenStrainShell7() { }

    virtual double give(Dof *dof, ValueModeType mode, double time);

    virtual bcType giveType() const { return DirichletBT; }

    /**
     * Initializes receiver according to object description stored in input record.
     * The input record contains two fields;
     * - gradient \#rows \#columns { d_11 d_12 ... ; d_21 ... } (required)
     * - cCoords \#columns x_1 x_2 ... (optional, default 0)
     * The prescribed tensor's columns must be equal to the size of the center coordinates.
     * The size of the center coordinates must be equal to the size of the coordinates in the applied nodes.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    /**
     * Constructs a coefficient matrix for all prescribed unknowns.
     * Helper routine for computational homogenization.
     * @todo Perhaps this routine should only give results for the dof it prescribes?
     * @param C Coefficient matrix to fill.
     */
    //void updateCoefficientMatrix(FloatMatrix &C);


    void evalCovarBaseVectorsAt(FloatMatrix &gcov, FloatArray &genEps, double zeta);
    void evalInitialCovarBaseVectorsAt(FloatMatrix &Gcov, FloatArray &genEps, double zeta);
    void setDeformationGradient(double zeta);
    void evaluateHigherOrderContribution(FloatArray &answer, double zeta, FloatArray &dx);

    /**
     * Set the center coordinate for the prescribed values to be set for.
     * @param x Center coordinate.
     */
    virtual void setCenterCoordinate(const FloatArray &x) { centerCoord = x; }
    /// Returns the center coordinate
    virtual FloatArray &giveCenterCoordinate() { return centerCoord; }

    virtual const char *giveClassName() const { return "PrescribedGenStrainShell7"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGenStrainShell7_Name; }

protected:
    double domainSize();
};
} // end namespace oofem

#endif // PrescribedGenStrainShell7_h
