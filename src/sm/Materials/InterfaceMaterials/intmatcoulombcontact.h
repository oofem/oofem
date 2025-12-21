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

#ifndef inmatcoulombcontact_h
#define inmatcoulombcontact_h

#include "material.h"
#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatCoulombContact
//@{
#define _IFT_IntMatCoulombContact_Name "intmatcoulombcontact"
#define _IFT_IntMatCoulombContact_kn "kn"
#define _IFT_IntMatCoulombContact_frictCoeff "frictcoeff"
#define _IFT_IntMatCoulombContact_stiffCoeff "stiffcoeff"
#define _IFT_IntMatCoulombContact_normalClearance "normalclearance"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to IntMatCoulombContact.
 * @author Milan Jirisek
 * @author Vit Smilauer
 * @author Jim Brouzoulis
 */
class IntMatCoulombContactStatus : public StructuralInterfaceMaterialStatus
{
protected:
    FloatArrayF<2> shearStressShift, tempShearStressShift;

public:
    IntMatCoulombContactStatus(GaussPoint *g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    const FloatArrayF<2> &giveShearStressShift() const { return shearStressShift; }
    void setTempShearStressShift(const FloatArrayF<2> &newShearStressShift) { tempShearStressShift = newShearStressShift; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * This class represents a "simple" interface material which is linear elastic in 
 * the normal direction. It is similar to a Coulomb contact model in the sense that 
 * the maximum shear stress is bounded by the normal stress times the coefficient 
 * of friction.
 * 
 * @author Milan Jirisek
 * @author Vit Smilauer
 * @author Jim Brouzoulis
 */
class IntMatCoulombContact : public StructuralInterfaceMaterial
{
protected:
    // Normal stiffness
    double kn = 0.;
    // Reduction factor of normal stiffness when material goes to tension
    double stiffCoeff = 0.;
    // Friction coefficient
    double frictCoeff = 0.;
    /// Normal distance which needs to be closed when interface element should act in compression (distance is 0 by default).
    double normalClearance = 0.;

public:
    IntMatCoulombContact( int n, Domain *d );

    const char *giveInputRecordName() const override { return _IFT_IntMatCoulombContact_Name; }
    const char *giveClassName() const override { return "IntMatCoulombContact"; }

    FloatArrayF<3> giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<3,3> give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;

    std::pair<double, FloatArrayF<2>> computeEngTraction(const FloatArrayF<2> &tempShearStressShift, double normalJump, const FloatArrayF<2> &shearJump );

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    bool hasAnalyticalTangentStiffness() const override { return true; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<IntMatCoulombContactStatus>(gp); }
};
} // end namespace oofem
#endif // simpleinterfacemat_h
