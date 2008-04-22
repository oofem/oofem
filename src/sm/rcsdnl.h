/* $Header: /home/cvs/bp/oofem/sm/src/rcsdnl.h,v 1.8 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *********************************************************************************************
//   *** CLASS NONLOCAL ROTATING SMEARED CRACK MODEL WITH TRANSITION TO SCALAR DAMAGE ************
//   *********************************************************************************************

#ifndef rcsdnl_h
#define rcsdnl_h

#include "rcsde.h"
#include "structuralnonlocalmaterialext.h"
#include "cltypes.h"

class GaussPoint;


class RCSDNLMaterialStatus : public RCSDEMaterialStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
    /*
     * This class implements associated Material Status to RCSDNLMaterial.
     *
     */

protected:

    FloatArray nonlocalStrainVector, tempNonlocalStrainVector, localStrainVectorForAverage;

public:

    RCSDNLMaterialStatus(int n, Domain *d, GaussPoint *g);
    ~RCSDNLMaterialStatus();

    void   printOutputAt(FILE *file, TimeStep *tStep);

    const FloatArray &giveNonlocalStrainVector()            { return nonlocalStrainVector; }
    const FloatArray &giveTempNonlocalStrainVector()        { return tempNonlocalStrainVector; }
    void   setTempNonlocalStrainVector(const FloatArray &ls) { tempNonlocalStrainVector = ls; }

    const FloatArray &giveLocalStrainVectorForAverage()     { return localStrainVectorForAverage; }
    void   setLocalStrainVectorForAverage(const FloatArray &ls) { localStrainVectorForAverage = ls; }

    // definition
    const char *giveClassName() const { return "RCSDNLMaterialStatus"; }
    classType             giveClassID() const { return RCSDNLMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    // saves current context(state) into stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Interface requesting service.
     * In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model status is derived both from class representing local status and from class
     * NonlocalMaterialStatusExtension or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning adress of component or using pointer conversion from receiver to base class
     * NonlocalMaterialStatusExtension. If no nonlocal extension exists, NULL pointer is returned.
     */
    virtual Interface *giveInterface(InterfaceType);
};



class RCSDNLMaterial : public RCSDEMaterial, public StructuralNonlocalMaterialExtensionInterface
{
    /*
     *
     * DESCRIPTION
     * This class implements a Nonlocal Rotating Crack Model with transition to scalar damage
     * for fracture in smeared fashion
     * ( only material stiffness modification is required, no changes in
     * mesh topology).
     * Model according to Milan Jirasek RC-SD model.
     *
     */

protected:
    /**
     * Nondimensional parameter controlling the transition between rc and sd model,
     * with respect to shear stifness degradation.
     */
    double SDTransitionCoeff2;
    /**
     * Interaction radius, related to the nonlocal characteristic length of material.
     */
    double R;
    /**
     * Strain at complete failure. For exponential law, ef is the strain at the intersection
     * of the horizontal axis with the tangent to the softening curve at peak stress.
     */
    double ef;

public:

    RCSDNLMaterial(int n, Domain *d);
    ~RCSDNLMaterial();

    // identification and auxiliary functions
    const char *giveClassName() const { return "RCSDNLMaterial"; }
    classType giveClassID()         const { return RCSDNLMaterialClass; }

    /** Interface requesting service */
    virtual Interface *giveInterface(InterfaceType);
    // contextIOResultType    saveContext (FILE* stream, void *obj = NULL);
    // contextIOResultType    restoreContext(FILE* stream, void *obj = NULL);
    IRResultType initializeFrom(InputRecord *ir);
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    //virtual     void   updateStatusForNewCrack( GaussPoint*, int, double);

    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    /**
     * Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     * @param strainVector total strain vector in given integration point.
     * @param gp integration point to update.
     * @param atTime solution step indicating time of update.
     */
    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime);

    /**
     * Computes the value of nonlocal weight function in given point.
     * @param src coordinates of source point.
     * @param coord coordinates of point, where nonlocal weight function is evaluated.
     * @return value of weight function.
     */
    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);
    /**
     * Determines, whether receiver has bounded weighting function (limited support)
     * @return true if weighting function bounded, zero otherwise
     */
    virtual int hasBoundedSupport() { return 1; }
    /**
     * Determines the width (radius) of limited support of weighting function
     */
    virtual void giveSupportRadius(double &radius) { radius = this->R; }


#ifdef __PARALLEL_MODE
    /**
     * Updates domain before nonloc average (using updateDomainBeforeNonlocAverage service)
     * to ensure, that the localStrainVectorForAverage variable is correctly updated and
     * pack this localStrainVectorForAverage into given buffer.
     * @see Material::packUnknowns for description.
     * @param buff communication buffer
     * @param stepN solution step
     * @param ip integration point
     */
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    /**
     * Unpack localStrainVectorForAverage value from given buffer.
     * @see Material::unpackAndUpdateUnknowns service.
     * @param buff communication buffer
     * @param stepN solution step.
     * @param ip integration point
     */
    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    /**
     * Estimates the necessary pack size to hold all packed data of receiver.
     */
    int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip);
#endif

    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RCSDNLMaterialStatus(1, RCSDEMaterial :: domain, gp); }

protected:
    double giveCharacteristicElementLenght(GaussPoint *, const FloatArray &) { return 1.0; }
    double giveMinCrackStrainsForFullyOpenCrack(GaussPoint *gp, int i);
    double computeStrength(GaussPoint *, double)  { return this->Ft; }

    ////
};


#endif // rcsdnl_h




