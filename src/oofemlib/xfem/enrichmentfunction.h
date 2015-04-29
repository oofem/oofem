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

#ifndef enrichmentfunction_h
#define enrichmentfunction_h

#include "intarray.h"
#include "femcmpnn.h"

#define _IFT_DiscontinuousFunction_Name "discontinuousfunction"
#define _IFT_HeavisideFunction_Name "heavisidefunction"
#define _IFT_RampFunction_Name "rampfunction"

namespace oofem {
class EnrichmentItem;
class EnrichmentDomain;
class BasicGeometry;
class GaussPoint;

/**
 * Abstract class representing global shape function
 * Base class declares abstract interface common to all enrichment functions.
 * Particularly, evaluateEnrFuncAt() and evaluateEnrFuncDerivAt()
 * evaluate the value and spatial derivatives at a given point.
 * @author chamrova
 * @author Jim Brouzoulis
 * @author Erik Svenning
 */
class OOFEM_EXPORT EnrichmentFunction : public FEMComponent
{
public:
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param aDomain Reference to domain.
     */
    EnrichmentFunction(int n, Domain *aDomain) : FEMComponent(n, aDomain) { }
    /// Destructor
    virtual ~EnrichmentFunction() { }

    // New interface
    virtual void evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const = 0;
    virtual void evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet) const = 0;

    /**
     * Returns the discontinuous jump in the enrichment function when the lvel set function
     * changes sign, e.g. 1.0 for Heaviside, 2.0 for Sign and 0.0 for abs enrichment.
     * Used for combination of cohesive zones and XFEM.
     */
    virtual void giveJump(std :: vector< double > &oJumps) const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "EnrichmentFunction"; }
    int giveNumberOfDofs() const { return numberOfDofs; }

protected:

    int numberOfDofs;
};

/** Class representing Sign EnrichmentFunction. */
class OOFEM_EXPORT DiscontinuousFunction : public EnrichmentFunction
{
public:
    DiscontinuousFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 1;
    }

    virtual void evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const;
    virtual void evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet) const;

    virtual void giveJump(std :: vector< double > &oJumps) const {
        oJumps.clear();
        oJumps.push_back(2.0);
    }

    virtual const char *giveClassName() const { return "DiscontinuousFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_DiscontinuousFunction_Name; }
};

/**
 * Class representing Heaviside EnrichmentFunction.
 *
 * @author Erik Svenning
 */
class OOFEM_EXPORT HeavisideFunction : public EnrichmentFunction
{
public:
    HeavisideFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 1;
    }

    virtual void evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const;
    virtual void evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet) const;

    virtual void giveJump(std :: vector< double > &oJumps) const {
        oJumps.clear();
        oJumps.push_back(1.0);
    }

    virtual const char *giveClassName() const { return "HeavisideFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_HeavisideFunction_Name; }
};

/** Class representing the four classical linear elastic branch functions. */
class OOFEM_EXPORT LinElBranchFunction
{
public:
    LinElBranchFunction() { }
    virtual ~LinElBranchFunction() { }

    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const double &iR, const double &iTheta) const;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const double &iR, const double &iTheta) const;

    virtual void giveJump(std :: vector< double > &oJumps) const;
    virtual void giveJump(std :: vector< double > &oJumps, const double &iRadius) const;

    std :: string errorInfo(const char *func) const { return std :: string( giveClassName() ) + func; }
    virtual const char *giveClassName() const { return "LinElBranchFunction"; }
};

/** Class representing a branch function for cohesive cracks. */
class OOFEM_EXPORT CohesiveBranchFunction
{
public:
    CohesiveBranchFunction():mExponent(0.5) { }
    virtual ~CohesiveBranchFunction() { }

    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const double &iR, const double &iTheta) const;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const double &iR, const double &iTheta) const;

    virtual void giveJump(std :: vector< double > &oJumps) const;
    virtual void giveJump(std :: vector< double > &oJumps, const double &iRadius) const;

    std :: string errorInfo(const char *func) const { return std :: string( giveClassName() ) + func; }
    virtual const char *giveClassName() const { return "CohesiveBranchFunction"; }

protected:
    double mExponent;
};

/** Class representing bimaterial interface. */
class OOFEM_EXPORT RampFunction : public EnrichmentFunction
{
public:

    RampFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 1;
    }

    virtual void evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const;
    virtual void evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet) const;

    virtual void giveJump(std :: vector< double > &oJumps) const {
        oJumps.clear();
        oJumps.push_back(0.0);
    }

    virtual const char *giveClassName() const { return "RampFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_RampFunction_Name; }
};
} // end namespace oofem
#endif  // enrichmentfunction_h
