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

#ifndef mat_cebfip90_h
#define mat_cebfip90_h

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for CebFipSlip90Material
//@{
#define _IFT_CebFipSlip90Material_Name "cebfipslip90"
#define _IFT_CebFipSlip90Material_tmax "tmax"
#define _IFT_CebFipSlip90Material_tres "tres"
#define _IFT_CebFipSlip90Material_s1 "s1"
#define _IFT_CebFipSlip90Material_s2 "s2"
#define _IFT_CebFipSlip90Material_s3 "s3"
//@}

namespace oofem {

/**
 * This class implements associated Material Status to IsoInterfaceDamageMaterial.
 * Stores scalar measure of the largest strain level reached.
 */
class CebFipSlip90MaterialStatus : public StructuralInterfaceMaterialStatus
{
protected:
    /// Scalar measure of the largest slip reached in material.
    double kappa = 0.;
    /// Non-equilibrated scalar of the largest slip displacement.
    double tempKappa = 0.;

public:
    /// Constructor.
    CebFipSlip90MaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() const { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() const { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }

    const char *giveClassName() const override { return "CebFipSlip90MaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relation damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class CebFipSlip90Material : public StructuralInterfaceMaterial
{
protected:
    /// Max force (stress).
    double tmax = 0.;
    /// Slip valu at begining of yield plateau.
    double s1 = 0.;
    /// Slip at end of plateau.
    double s2 = 0.;
    /// Slip when residual force/stress activated.
    double s3 = 0.;
    /// Residual force/stress.
    double tres = 0.;
    /// Alpha coeff.
    double alpha = 0.;

public:
    /// Constructor
    CebFipSlip90Material(int n, Domain * d);

    const char *giveInputRecordName() const override { return _IFT_CebFipSlip90Material_Name; }
    const char *giveClassName() const override { return "CebFipSlip90Material"; }

    double giveEngTraction_1d(double jump, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<1,1> give1dStiffnessMatrix_Eng(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    bool hasAnalyticalTangentStiffness() const override { return true; }

    /**
     * Computes the value of bond force/stress, based on given value of slip value.
     * @param kappa Slip value.
     */
    double computeBondForce(double kappa) const;
    /**
     * Computes the value of bond force/stress stiffness, based on given value of slip value.
     * @param kappa Slip value.
     */
    double computeBondForceStiffness(double kappa) const;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<CebFipSlip90MaterialStatus>(gp); }
};
} // end namespace oofem
#endif // mat_cebfip90_h
