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

#ifndef mixedpressurematerialextensioninterface_h
#define mixedpressurematerialextensioninterface_h

#include "interface.h"
#include "matresponsemode.h"
#include "domain.h"

///@name micromorphicmaterialextensioninterface
//@{

//@}

namespace oofem {
class FloatMatrix;
class FloatArray;
class GaussPoint;
class TimeStep;



/**
 * Material interface for gradient material models.
 */
class MixedPressureMaterialExtensionInterface : public Interface
{
protected:
    Domain *dom = nullptr;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    MixedPressureMaterialExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~MixedPressureMaterialExtensionInterface() { }


    virtual void giveDeviatoric3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                         MatResponseMode,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep) const
    { OOFEM_ERROR("not implemented "); }


    virtual void giveDeviatoricPlaneStrainStiffMtrx(FloatMatrix &answer,
                                                    MatResponseMode, GaussPoint *gp,
                                                    TimeStep *tStep) const
    { OOFEM_ERROR("not implemented "); }


    virtual void giveDeviatoricConstitutiveMatrix(FloatMatrix & answer,
                                                  MatResponseMode, GaussPoint * gp,
                                                  TimeStep * tStep) const;


    virtual void giveInverseOfBulkModulus(double &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;


    void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep) const;

    virtual void giveRealStressVectorUP_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep) const = 0;
    virtual void giveRealStressVectorUP_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep) const;
    //    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep) = 0;

    virtual void giveFiniteStrainGeneralizedStressVectors(FloatArray &sigma, GaussPoint *gp, const FloatArray &devF, double pressure, TimeStep *tStep) {; }
};
}
#endif
