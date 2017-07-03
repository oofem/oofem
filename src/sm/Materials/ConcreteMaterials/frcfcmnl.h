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

#ifndef frcfcmnl_h
#define frcfcmnl_h


#include "frcfcm.h"
#include "Materials/structuralnonlocalmaterialext.h"

///@name Input fields for frcfcmnl
//@{
#define _IFT_FRCFCMNL_Name "frcfcmnl"
#define _IFT_FRCFCMNL_participAngle "alpha"
//@}

namespace oofem {
class GaussPoint;
/**
 * This class implements a FRCFCMNL material in a finite element problem.
 *
 * It extends frcfcm model to behave as nonlocal in terms of crack spacing
 *
 */

class FRCFCMNLStatus : public FRCFCMStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:

    /// bulk stress in fibers - evaluated from crack opening
    FloatArray fiberStressLoc;
    /// Non-equilibrated stress (bulk) in fibers.
    FloatArray tempFiberStressLoc;

    /// bulk stress in fibers - evaluated from crack opening
    FloatArray fiberStressNL;
    /// Non-equilibrated stress (bulk) in fibers.
    FloatArray tempFiberStressNL;


public:
    FRCFCMNLStatus(int n, Domain *d, GaussPoint *g);
    virtual ~FRCFCMNLStatus();

    /// LOCAL FIBER STRESSES (from crack opening)
    double giveFiberStressLoc(int icrack) { return fiberStressLoc.at(icrack); }
    double giveTempFiberStressLoc(int icrack) { return tempFiberStressLoc.at(icrack); }
    void setTempFiberStressLoc(int icrack, double newFiberStressLoc) { tempFiberStressLoc.at(icrack) = newFiberStressLoc; }

    /// NON-LOCAL FIBER STRESSES (from surrounding cracks)
    double giveFiberStressNL(int icrack) { return fiberStressNL.at(icrack); }
    double giveTempFiberStressNL(int icrack) { return tempFiberStressNL.at(icrack); }
    void setTempFiberStressNL(int icrack, double newFiberStressNL) { tempFiberStressNL.at(icrack) = newFiberStressNL; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "FRCFCMNLStatus"; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);



    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);
};


class FRCFCMNL : public FRCFCM, public StructuralNonlocalMaterialExtensionInterface
{
public:
    FRCFCMNL(int n, Domain *d);
    virtual ~FRCFCMNL() {}

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "FRCFCMNL"; }
    virtual const char *giveInputRecordName() const { return _IFT_FRCFCMNL_Name; }

    // identification and auxiliary functions
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);


    virtual void giveMaterialStiffnessMatrix(FloatMatrix & answer, MatResponseMode,
                                             GaussPoint * gp,
                                             TimeStep * tStep);

    virtual Interface *giveInterface(InterfaceType it);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new FRCFCMNLStatus(1, FRCFCM :: domain, gp); }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    //nothing to update here, is it?
    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) {}

    virtual bool isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress);


    double computeDebondedLength(double delta);

    /// compute the the difference in fiber stress in the target (local stress) and receiver (nonlocal stress)
    double computeDecreaseInFibreStress(double distance, double delta, double debondedLength);

    /// computes cetroid of a finite element - works only for linear 3 and 4-node elements
    void computeElementCentroid(FloatArray &answer, GaussPoint *gp);

    /// checks if a element center of homeGP is in projection of element containing nearGP
    bool isInElementProjection(GaussPoint *homeGp, GaussPoint *nearGp, int iNlCrack);

    /// computes nonlocal stress in fibers in cracked GP
    virtual double computeNonlocalStressInFibers(const FloatArray &crackVector, GaussPoint *gp, TimeStep *tStep);

    /// computes nonlocal stress in fibers in uncracked GP
    virtual double computeNonlocalStressInFibersInUncracked(GaussPoint *gp, TimeStep *tStep);

    /// computes an angle between two vectors
    double computeAngleBetweenVectors(const FloatArray &vec1, const FloatArray &vec2);


protected:
    /// participation angle. The target gauss point must fall into this angle to contribute to the nonlocal stress
    double participAngle;
};
} // end namespace oofem
#endif // frcfcmnl_h
