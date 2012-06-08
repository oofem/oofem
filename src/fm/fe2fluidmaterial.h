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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef fe2sinteringmaterial_h
#define fe2sinteringmaterial_h

#include "fluiddynamicmaterial.h"
#include "matstatus.h"
#include "mixedgradientpressurebc.h"

namespace oofem {

class StokesFlow;

/**
 * Class representing material status for the subscale fluid, i.e an Representative Volume Element (RVE).
 * @author Mikael Öhman
 */
class FE2FluidMaterialStatus : public FluidDynamicMaterialStatus
{
protected:
    /// The subscale flow
    StokesFlow *rve;
    /// Boundary condition in RVE that performs the computational homogenization.
    MixedGradientPressureBC *bc;

    FloatMatrix Ed;
    FloatArray Cd;
    FloatArray Ep;
    double Cp;

    double voffraction;

public:
    /**
     * Creates new material status.
     * @param n Material status number.
     * @param d Domain that status belongs to.
     * @param gp Gauss point that the status belongs to.
     * @param inputfile The input file describing the micro problem.
     * @param porosity Initial porosity.
     */
    FE2FluidMaterialStatus(int n, Domain *d, GaussPoint *gp, const std::string &inputfile);
    /// Destructor
    virtual ~FE2FluidMaterialStatus();

    StokesFlow *giveRVE() { return this->rve; }
    MixedGradientPressureBC *giveBC() { return this->bc; }

    double giveVOFFraction() { return this->voffraction; }

    /// Creates/Initiates the RVE problem.
    virtual bool createRVE(int n, GaussPoint *gp, const std::string &inputfile);

    /// Copies time step data to RVE.
    virtual void setTimeStep(TimeStep *tStep);

    double computeSize();

    FloatMatrix &giveDeviatoricTangent() { return Ed; }
    FloatArray &giveDeviatoricPressureTangent() { return Ep; };
    FloatArray &giveVolumetricDeviatoricTangent() { return Cd; };
    double &giveVolumetricPressureTangent() { return Cp; };

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "FE2FluidMaterialStatus"; }
    virtual classType giveClassID() const { return FE2FluidMaterialStatusClass; }
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
    std::string inputfile;
    static int n;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    FE2FluidMaterial(int n, Domain *d) : FluidDynamicMaterial(n, d) { }
    /// Destructor.
    virtual ~FE2FluidMaterial() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual int checkConsistency();
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void computeDeviatoricStressVector(FloatArray &stress_dev, double &r_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep);
    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveDeviatoricPressureStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveVolumetricDeviatoricStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveVolumetricPressureStiffness(double &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual const char *giveClassName() const { return "FE2FluidMaterial"; }
    virtual classType giveClassID() const { return FE2FluidMaterialClass; }
};

} // end namespace oofem
#endif // rvesinteringmaterial_h
