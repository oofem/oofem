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

#ifndef latticebondplasticity_h
#define latticebondplasticity_h

#include "latticelinearelastic.h"
#include "latticematstatus.h"

///@name Input fields for LatticePlasticityDamage
//@{
#define _IFT_LatticeBondPlasticity_Name "latticebondplast"
#define _IFT_LatticeBondPlasticity_tol "tol"
#define _IFT_LatticeBondPlasticity_iter "iter"
#define _IFT_LatticeBondPlasticity_sub "sub"
#define _IFT_LatticeBondPlasticity_fc "fc"
#define _IFT_LatticeBondPlasticity_angle1 "angle1"
//#define _IFT_LatticeBondPlasticity_angle2 "angle2"
//#define _IFT_LatticeBondPlasticity_flow "flow"
#define _IFT_LatticeBondPlasticity_ef "ef"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to LatticeBondPlasticity.
 */
class LatticeBondPlasticityStatus : public LatticeMaterialStatus
{
protected:
    /* FloatArray plasticStrain; */
    /* FloatArray tempPlasticStrain; */

    double kappaP;
    double tempKappaP;

    int surfaceValue;

public:

    /// Constructor
    LatticeBondPlasticityStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~LatticeBondPlasticityStatus() {}


    /// Returns the last equilibrated scalar measure of the largest strain level
    double giveKappaP() { return kappaP; }
    /// Returns the temp. scalar measure of the largest strain level
    double giveTempKappaP() { return tempKappaP; }
    void  letTempPlasticStrainBe(const FloatArray &v)
    { tempPlasticStrain = v; }
    /// Sets the temp scalar measure of the largest strain level to given value
    void setTempKappaP(double newKappa) { tempKappaP = newKappa; }
    void setSurfaceValue(int val)
    {
        surfaceValue = val;
    }

    int giveSurfaceValue()
    {
        return surfaceValue;
    }


    virtual void   printOutputAt(FILE *file, TimeStep *tStep);

    const char *giveClassName() const { return "LatticeBondPlasticityStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);


    void  saveContext(DataStream &stream, ContextMode mode);

    void  restoreContext(DataStream &stream, ContextMode mode);
};

class LatticeBondPlasticity : public LatticeLinearElastic
{
protected:

    enum LatticeBondPlasticity_SurfaceType { ST_Vertex, ST_Shear, ST_Compression };
    LatticeBondPlasticity_SurfaceType surfaceType;

    enum LatticeBondPlasticity_ReturnResult { RR_NotConverged, RR_Converged, RR_Elastic};
    LatticeBondPlasticity_ReturnResult returnResult;

    double initialYieldStress;

    /// compressive strength
    double fc;
    /// frictional angle of the yield surface
    double frictionAngleOne;
    double frictionAngleTwo;

    double hardeningLimit;

    /// frictional angle of the plastic potential
    double flowAngle;

    /// determines the softening -> corresponds to crack opening (not strain) when tension stress vanishes
    double wf;

    int oldApproachFlag;

    //pi constant
    double myPi;

    /// yield tolerance
    double yieldTol;

    /// maximum number of iterations for stress return
    int newtonIter;

    /// maximum number of subincrements
    int numberOfSubIncrements;

    //parameter in hardening law
    double ef;
    
public:

    /// Constructor
    LatticeBondPlasticity(int n, Domain *d);
    /// Destructor
    ~LatticeBondPlasticity();

    virtual const char *giveInputRecordName() const { return _IFT_LatticeBondPlasticity_Name; }
    const char *giveClassName() const { return "LatticeBondPlasticity"; }

    double computeHardening(double kappa);
    
    double computeDHardeningDKappa(double kappa);

    double computeParamA(const double kappa);
    double computeDParamADKappa(const double kappa);

    double computeShift(const double kappa);
    double computeDShiftDKappa(const double kappa);

      
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }



    bool checkForVertexCase(FloatArray &stress, GaussPoint *gp);

    void performVertexReturn(FloatArray &stress, GaussPoint *gp);

    double computeTransition(const double kappa, GaussPoint *gp);

    int hasMaterialModeCapability(MaterialMode mode);


    /* Compute the A matrix used for closest point return
     */
    void computeBMatrix(FloatMatrix &answer,
                        const FloatArray &sigma,
                        const double deltaLambda);
    void computeAMatrix(FloatMatrix &answer,
                        const FloatArray &sigma,
                        const double tempKappa,
                        const double deltaLambda,
                        int transitionFlag,
                        GaussPoint *gp);


    void computeJacobian(FloatMatrix &answer,
                         const FloatArray &sigma,
                         const double tempKappa,
                         const double deltaLambda,
                         int transitionFlag,
                         GaussPoint *gp);

    int computeInverseOfJacobian(FloatMatrix &answer,
                                 const FloatMatrix &src);

    void computeFVector(FloatArray &answer,
                        const FloatArray &sigma,
                        const double kappa,
                        int transitionFlag,
                        GaussPoint *gp);
    void computeMVector(FloatArray &answer,
                        const FloatArray &sigma,
                        const double kappa,
                        int transitionFlag,
                        GaussPoint *gp);
    void computeDMMatrix(FloatMatrix &answer,
                         const FloatArray &sigma,
                         const double deltaLambda,
                         int transitionFlag,
                         GaussPoint *gp);


    void giveRealStressVector(FloatArray &answer, GaussPoint *,
                              const FloatArray &, TimeStep *);


    void performPlasticityReturn(FloatArray &answer,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain, TimeStep *);

    double performRegularReturn(FloatArray &stress,
                                double yieldValue,
                                int transitionFlag,
                                GaussPoint *gp);

    double computeYieldValue(const FloatArray &sigma,
                             const double tempKappa,
                             int transitionFlag,
                             GaussPoint *gp);


    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void giveReducedStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

protected:
};
} // end namespace oofem
#endif
