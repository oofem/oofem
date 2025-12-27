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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef isoasymm1d_h
#define isoasymm1d_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"


///@name Input fields 
//@{
#define _IFT_IsotropicAsymmetric1DMaterial_Name "isoasymm1d"
#define _IFT_IsotropicAsymmetric1DMaterial_ec "ec"
#define _IFT_IsotropicAsymmetric1DMaterial_et "et"
#define _IFT_IsotropicAsymmetric1DMaterial_efc "efc"
#define _IFT_IsotropicAsymmetric1DMaterial_eft "eft"
#define _IFT_IsotropicAsymmetric1DMaterial_talpha "talpha"
#define _IFT_IsotropicAsymmetric1DMaterial_m "m"

//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements an isotropic asymmetric (different stiffness in tension and compression) linear elastic material in a finite element problem.
 *
 * Tasks:
 * - Returning standard material stiffness matrix for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint.
 * - methods give2dPlaneStressMtrx, givePlaneStrainMtrx, give1dStressMtrx are
 *   introduced since form of this matrices is well known, and for
 *   faster response mainly in linear elastic problems.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class IsotropicAsymmetric1DMaterial : public StructuralMaterial
{
protected:
    /// Young's modulus in tension.
    double Et = 0;
    /// Young's modulus in compression
    double Ec = 0;
    // failure strain in compression 
    double efc = 1.0; // positive value means efc not set, no failure in compression
    // failure strain in tension
    double eft = -1.0; // negative values means eft not set, no failure in tension
    /// Alpha
    double a = 0;
    /// Regularization parameter
    double m = 15;

public:
    /**
     * Creates a new IsotropicAsymmetric1DMaterial class instance
     * with given number belonging to domain d.
     * @param n material model number in domain
     * @param d domain which receiver belongs to
     */
    IsotropicAsymmetric1DMaterial(int n, Domain *d);
    /**
     * Creates a new IsotropicAsymmetric1DMaterial class instance
     * with given number belonging to domain d.
     * @param n Material model number in domain.
     * @param d Domain which receiver belongs to.
     * @param E Young modulus.
     * @param nu Poisson ratio.
     */
    IsotropicAsymmetric1DMaterial(int n, Domain *d, double Et, double Ec, double efc, double eft);

    const char *giveClassName() const override { return "IsotropicAsymmetric1DMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_IsotropicAsymmetric1DMaterial_Name; }
    bool hasMaterialModeCapability(MaterialMode mode) const override;

    /**
     * Initializes receiver according to object description stored in input record.
     * The E modulus (keyword "E"), Poisson ratio ("nu") and coefficient of thermal dilatation
     * alpha ("talpha") are read. The parent class instanciateFrom method is called.
     */
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double give(int aProperty, GaussPoint *gp) const override;

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode, GaussPoint * gp,
                                            TimeStep * tStep) const override;
    virtual FloatArrayF< 1 >giveRealStressVector_1d(const FloatArrayF< 1 > &reducedE, GaussPoint *gp, TimeStep *tStep) const override;
    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<StructuralMaterialStatus>(gp); }

};
} // end namespace oofem
#endif // isoasymm1d_h
