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

#ifndef graddpmaterialextensioninterface_h
#define graddpmaterialextensioninterface_h

#include "interface.h"
#include "matresponsemode.h"

///@name graddpdmaterialextensioninterface
//@{
#define _IFT_GradDpMaterialExtensionInterface_averagingtype "averagingType"
#define _IFT_GradDpMaterialExtensionInterface_cl "cl"
#define _IFT_GradDpMaterialExtensionInterface_beta "beta"
#define _IFT_GradDpMaterialExtensionInterface_zeta "zeta"
//@}

namespace oofem {
class FloatMatrix;
class GaussPoint;
class TimeStep;



/**
 * Material interface for gradient material models.
 */
class GradDpMaterialExtensionInterface : public Interface
{
protected:
    Domain *dom;

    /**
     * Initial(user defined) characteristic length of the nonlocal model
     * (its interpretation depends on the weight function)
     * Is different to cl when a Stress-based or a Distance-based
     * nonlocal variation is applied
     */
    double cl0;


    /**
     * Parameter which defines the averaging type
     * When averType is equal to zereo classical approach is used
     * for averType equal to one, distance-based apprach is used
     * and avetType equal to two corresponds to stress-based averaging
     */
    int averType;

    /**
     * Parameter which multiplied with the interaction radius cl0
     * gives its minimum allowed value. It is used when a Stress-based
     * or a Distance-based nonlocal variation is applied
     */
    double beta;
    /**
     * Parameter used when Distance-based nonlocal variation is applied
     * When it is multiplied with the interaction radius cl gives the maxinmum
     * distance of the Gauss Point from the boundary. If the Gauss Point's distance
     * from the boundary is larger than this value the interaction radius cl is set
     * to cl0
     */
    double zeta;

    /**
     * Characteristic length of the nonlocal model
     * (its interpretation depends on the type of weight function).
     */
    double cl;

    /**
     * Provides the distance based interaction radius
     * This function is called when averType is set to 1.
     * The function loops over all user defined nonlocal boundaries
     * to find minimum distance from the GP. Then calculates interaction radius
     * @param gpCoords The Gauss points' coordinates, whose interaction radius is calculated based on the distance-based averaging approach.
     */
    void giveDistanceBasedCharacteristicLength(const FloatArray &gpCoords);
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    GradDpMaterialExtensionInterface(Domain * d);
    /// Destructor.
    virtual ~GradDpMaterialExtensionInterface() { }
    /// Left upper block
    virtual void givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Left lower block
    virtual void givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Right upper block
    virtual void givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Right lower block
    virtual void givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// Stress-based averaging
    virtual void givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /// gradient - based giveRealStressVector
    virtual void giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivningVariable, TimeStep *tStep) { OOFEM_ERROR("not implemented") }
    virtual void giveFirstPKStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivningVariable, TimeStep *tStep) { OOFEM_ERROR("not implemented") }
    virtual void giveCauchyStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalDamageDrivningVariable, TimeStep *tStep) { OOFEM_ERROR("not implemented") }

    virtual IRResultType initializeFrom(InputRecord *ir);
    int giveAveragingType() { return averType; }
};

class GradDpMaterialStatusExtensionInterface : public Interface
{
protected:
    double nonlocalCumulatedStrain;
public:
    virtual double giveNonlocalCumulatedStrain() = 0;
    virtual void setNonlocalCumulatedStrain(double nonlocalCumulatedStrain) = 0;
};
}
#endif
