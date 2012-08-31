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

#ifndef hydratingconcretemat_h
#define hydratingconcretemat_h

#include "isoheatmat.h"
#include "hydratingisoheatmat.h"

namespace oofem {
/**
 * This class implements various phenomenological and affinity hydration models. No coupling with relative humidity
 * is considered. Heat capacity and thermal conductivity can be set constant or concrete may be treated as a 5-component
 * evolving material.
 */
class HydratingConcreteMat : public IsotropicHeatTransferMaterial
{
public:
    HydratingConcreteMat(int n, Domain *d);
    virtual ~HydratingConcreteMat();

    /// Return true if hydration heat source is present.
    virtual int hasInternalSource() { return 1; };
    virtual void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode);
    virtual void updateInternalState(const FloatArray &state, GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "HydratingConcreteMat"; }
    virtual classType giveClassID() const { return HydratingConcreteMatClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    // post-processing
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
    virtual double giveConcreteConductivity(GaussPoint *gp);
    virtual double giveConcreteCapacity(GaussPoint *gp);
    virtual double giveConcreteDensity(GaussPoint *gp);
    /// Type of hydration model, e.g. exponential curve, Cervera's model.
    int hydrationModelType;
    double maxModelIntegrationTime;
    /// Minimum number of integration steps for hydration model within a given timeStep.
    double minModelTimeStepIntegrations;
    /// Potential heat of hydration, for ordinary Portland cement approximately 500 J/g.
    double Qpot;
    /// Mass of cement in kg per 1m3 of concrete.
    double massCement;
    /// Activation energy of concrete (default 38400 J/mol/K).
    double activationEnergy;
    /// Reference temperature for hydration model.
    double referenceTemperature;
    /**
     * Parameters for exponential affinity hydration model summarized in A.K. Schindler and K.J. Folliard:
     * Heat of Hydration Models for Cementitious Materials, ACI Materials Journal, 2005.
     */
    double tau, beta;

    /**
     * Parameters for affinity hydration model inspired by Cervera et al.
     * Journal of Engineering Mechanics ASCE, 125(9), 1018-1027, 1999.
     */
    double B1, B2, eta, DoHInf;

protected:
    /// Use different methods to evaluate material conductivity, capacity, or density
    int conductivityType, capacityType, densityType;
    /// Degree of reinforcement, if defined, reinforcement effect for conductivity and capacity is accounted for. Isotropic case.
    double reinforcementDegree;
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};

/**
 * HydratingConcreteMatStatus stores degree of hydration in each integration point
 */
class HydratingConcreteMatStatus : public TransportMaterialStatus
{
public:
    HydratingConcreteMatStatus(int n, Domain *d, GaussPoint *g);
    virtual ~HydratingConcreteMatStatus();
    double power;
    double GivePower(TimeStep *atTime);
    /// Returns actual degree of hydration at last known equilibrium.
    virtual double giveDoHActual(void);
    virtual void updateYourself(TimeStep *atTime);
    virtual void printOutputAt(FILE *file, TimeStep *atTime);
    double lastIntrinsicTime;
protected:
    double lastEquivalentTime, equivalentTime, degreeOfHydration, lastDegreeOfHydration;
    double scaleTemperature(void);
    /// Return affinity scaled to 25C.
    double affinity25(double alpha);
};
} // end namespace oofem
#endif // hydratingconcretemat_h
