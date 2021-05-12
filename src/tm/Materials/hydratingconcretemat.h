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

#ifndef hydratingconcretemat_h
#define hydratingconcretemat_h

#include "tm/Materials/isoheatmat.h"
#include "tm/Materials/hydratingisoheatmat.h"

///@name Input fields for HydratingConcreteMat
//@{
#define _IFT_HydratingConcreteMat_Name "hydratingconcretemat"
#define _IFT_HydratingConcreteMat_referenceTemperature "referencetemperature"
#define _IFT_HydratingConcreteMat_castAt "castat"
#define _IFT_HydratingConcreteMat_hydrationModelType "hydrationmodeltype"
#define _IFT_HydratingConcreteMat_maxModelIntegrationTime "maxmodelintegrationtime"
#define _IFT_HydratingConcreteMat_minModelTimeStepIntegrations "minmodeltimestepintegrations"
#define _IFT_HydratingConcreteMat_conductivitytype "conductivitytype"
#define _IFT_HydratingConcreteMat_capacitytype "capacitytype"
#define _IFT_HydratingConcreteMat_densitytype "densitytype"
#define _IFT_HydratingConcreteMat_activationEnergy "activationenergy"
#define _IFT_HydratingConcreteMat_massCement "masscement"
#define _IFT_HydratingConcreteMat_reinforcementDegree "reinforcementdegree"
#define _IFT_HydratingConcreteMat_tau "tau"
#define _IFT_HydratingConcreteMat_beta "beta"
#define _IFT_HydratingConcreteMat_B1 "b1"
#define _IFT_HydratingConcreteMat_B2 "b2"
#define _IFT_HydratingConcreteMat_eta "eta"
#define _IFT_HydratingConcreteMat_DoHInf "dohinf"
#define _IFT_HydratingConcreteMat_DoH1 "doh1"
#define _IFT_HydratingConcreteMat_P1 "p1"
#define _IFT_HydratingConcreteMat_qpot "qpot"
#define _IFT_HydratingConcreteMat_wc "w/c"
#define _IFT_HydratingConcreteMat_ac "a/c"
#define _IFT_HydratingConcreteMat_rhoCem "rhocem"
#define _IFT_HydratingConcreteMat_rhoAgg "rhoagg"
#define _IFT_HydratingConcreteMat_Blaine "blaine"
#define _IFT_HydratingConcreteMat_alphaSet0 "alphaset0"
#define _IFT_HydratingConcreteMat_timeSet "timeset"
#define _IFT_HydratingConcreteMat_alphaCrit0 "alphacrit0"
#define _IFT_HydratingConcreteMat_B0 "b0"
#define _IFT_HydratingConcreteMat_timeToSeconds "timetoseconds"


//@}

namespace oofem {
/**
 * This class implements various phenomenological and affinity hydration models. No coupling with relative humidity
 * is considered. Heat capacity and thermal conductivity can be set constant or concrete may be treated as a 5-component
 * evolving material.
 */
class HydratingConcreteMat : public IsotropicHeatTransferMaterial
{
public:
    HydratingConcreteMat(int n, Domain * d);

    bool hasInternalSource() const override { return true; }
    bool hasCastingTimeSupport() const override { return true; }
    void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    double giveCharacteristicValue(MatResponseMode mode,
                                   GaussPoint *gp,
                                   TimeStep *tStep) const override;

    const char *giveClassName() const override { return "HydratingConcreteMat"; }

    void initializeFrom(InputRecord &ir) override;

    // post-processing
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    double giveIsotropicConductivity(GaussPoint *gp, TimeStep *tStep) const override;
    virtual double giveConcreteCapacity(GaussPoint *gp, TimeStep *tStep) const;
    virtual double giveConcreteDensity(GaussPoint *gp, TimeStep *tStep) const;

    /// Type of hydration model, e.g. exponential curve, Cervera's model.
    int hydrationModelType = 0;
    double maxModelIntegrationTime = 0.;
    /// Minimum number of integration steps for hydration model within a given timeStep.
    double minModelTimeStepIntegrations = 0.;
    /// Potential heat of hydration, for ordinary Portland cement approximately 500 J/g.
    double Qpot = 0.;
    /// Mass of cement in kg per 1m3 of concrete.
    double massCement = 0.;
    /// Activation energy of concrete (default 38400 J/mol/K).
    double activationEnergy = 0.;
    /// Reference temperature for hydration model.
    double referenceTemperature = 0.;
    /**
     * Parameters for exponential affinity hydration model summarized in A.K. Schindler and K.J. Folliard:
     * Heat of Hydration Models for Cementitious Materials, ACI Materials Journal, 2005.
     */
    double tau = 0., beta = 0.;

    /**
     * Parameters for affinity hydration model inspired by Cervera et al.
     * Journal of Engineering Mechanics ASCE, 125(9), 1018-1027, 1999.
     */
    double B1 = 0., B2 = 0., eta = 0., DoHInf = 0.;
    ///Optional extension to slag-rich, high-blended cements

    double DoH1 = 0., P1=0.;
    /**
     * Parameters for hydration model Saeed Rahimi-Aghdam, Zdeněk P. Bažant, Gianluca Cusatis: Extended Microprestress-Solidification Theory (XMPS) for Long-Term Creep and Diffusion Size Effect in Concrete at Variable Environment, JEM-ASCE, 2019. Appendix A.
     */
    ///Water/cement ratio and aggregate/cement ratio  
    double wc, ac;
    ///Density of cement and aggregates (weighted average from fine and coarse aggregates
    double rhoCem, rhoAgg;
    ///Initial volume fraction of cement and water
    double Vc0, Vw0;
    ///Volume fractions at setting time
    double VCemSet, VCHSet, VGelSet;
    ///Average cement particle radius (m)
    double a0;
    ///Number of cement particles in a unit volume
    double ng;
    ///Degree of hydration for setting time
    double alphaSet;
    ///Time at setting
    double timeSet;
    ///Radius of cement particle at setting time
    double aSet;
    ///Radius of gel barrier at setting time
    double zSet;
    ///Degree of hydration at which the gel barrier will be completed
    double alphaCrit;
    ///Basic diffusivity (about 1.1e-11 m2/day for Portland cements)
    double B0;
    ///timeToSeconds, =1 when running in seconds (default), =86400 when time in days
    double timeToSeconds;

protected:
    double GivePower(TimeStep *tStep, GaussPoint *gp, ValueModeType mode) const;
    double scaleTemperature(GaussPoint *gp) const;
    /// Return affinity scaled to 25C.
    double affinity25(double alpha) const;

    /// Use different methods to evaluate material conductivity, capacity, or density
    int conductivityType, capacityType, densityType;
    /// Degree of reinforcement, if defined, reinforcement effect for conductivity and capacity is accounted for. Isotropic case.
    double reinforcementDegree;
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};

/**
 * HydratingConcreteMatStatus stores degree of hydration in each integration point
 */
class HydratingConcreteMatStatus : public TransportMaterialStatus
{
public:
    HydratingConcreteMatStatus(GaussPoint * g);
    /// Returns actual degree of hydration at last known equilibrium.
    double giveDoHActual() const;
    void updateYourself(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    double power = 0.;
    double lastEvalTime = -1.e20;
    double lastEquivalentTime = 0., equivalentTime = 0., degreeOfHydration = 0., lastDegreeOfHydration = 0.;
    /// Radius of the equivalent contact-free C-S-H shells
    double zShell = 0., lastZShell = 0.;
    // Radius of cement particle
    double aCement = 0., lastACement = 0.;
    // Volume fractions of cement, gel, CH;
    double VCem = 0., lastVCem = 0., VGel = 0., lastVGel = 0., VCH = 0., lastVCH = 0.;
};
} // end namespace oofem
#endif // hydratingconcretemat_h
