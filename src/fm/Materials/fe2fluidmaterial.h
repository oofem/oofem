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

#ifndef fe2sinteringmaterial_h
#define fe2sinteringmaterial_h

#include "fm/Materials/fluiddynamicmaterial.h"
#include "matstatus.h"
#include "mixedgradientpressurebc.h"
#include "floatmatrix.h"

#include <memory>

///@name Input fields for FE^2 fluid material
//@{
#define _IFT_FE2FluidMaterial_Name "fe2fluidmaterial"
#define _IFT_FE2FluidMaterial_fileName "inputfile"
//@}

namespace oofem {
class EngngModel;

/**
 * Class representing material status for the subscale fluid, i.e an Representative Volume Element (RVE).
 * @author Mikael Öhman
 */
class FE2FluidMaterialStatus : public FluidDynamicMaterialStatus
{
protected:
    /// The subscale flow
    std :: unique_ptr< EngngModel > rve;
    /// Boundary condition in RVE that performs the computational homogenization.
    MixedGradientPressureBC *bc;

    FloatMatrix Ed;
    FloatArray Cd;
    FloatArray Ep;
    double Cp;

    double pressure;
    double voffraction;

    bool oldTangents;

public:
    /**
     * Creates new material status.
     * @param n Material status number.
     * @param rank Parallel rank (ignored if negative).
     * @param gp Gauss point that the status belongs to.
     * @param inputfile The input file describing the micro problem.
     */
    FE2FluidMaterialStatus(int n, int rank, GaussPoint * gp, const std :: string & inputfile);

    EngngModel *giveRVE() { return this->rve.get(); }
    MixedGradientPressureBC *giveBC() { return this->bc; }

    void markOldTangents();
    void computeTangents(TimeStep *tStep);

    double giveVOFFraction() { return this->voffraction; }

    /// Creates/Initiates the RVE problem.
    bool createRVE(int n, int rank, GaussPoint *gp, const std :: string &inputfile);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    FloatMatrix &giveDeviatoricTangent() { return Ed; }
    FloatArray &giveDeviatoricPressureTangent() { return Ep; }
    FloatArray &giveVolumetricDeviatoricTangent() { return Cd; }
    double &giveVolumetricPressureTangent() { return Cp; }

    double givePressure() { return this->pressure; }
    void letPressureBe(double val) { this->pressure = val; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "FE2FluidMaterialStatus"; }
};


/**
 * Multiscale constitutive model for subscale flow problems, typically sintering.
 *
 * The material uses the MixedGradientPressureBC to perform computational homogenization.
 * The requirement for the supplied subscale flow problem is:
 * - It must have a MixedGradientPressureBC. It should be the only Dirichlet boundary condition.
 * - It must have boundary elements along the outer boundary such that the total volume (including internal pores) can be computed.
 *
 * @author Mikael Öhman
 */
class FE2FluidMaterial : public FluidDynamicMaterial
{
private:
    std :: string inputfile;
    static int n;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    FE2FluidMaterial(int n, Domain * d) : FluidDynamicMaterial(n, d) { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    int checkConsistency() override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    std::pair<FloatArrayF<6>, double> computeDeviatoricStress3D(const FloatArrayF<6> &eps, double pressure, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF<6> computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<6,6> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    Tangents<6> computeTangents3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    const char *giveClassName() const override { return "FE2FluidMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_FE2FluidMaterial_Name; }
};
} // end namespace oofem
#endif // rvesinteringmaterial_h
