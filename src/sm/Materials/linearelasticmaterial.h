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

#ifndef linearelasticmaterial_h
#define linearelasticmaterial_h

#include "../sm/Materials/structuralmaterial.h"

namespace oofem {
/**
 * This class is a abstract base class for all linear elastic material models
 * in a finite element problem. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 * Efficient implementation of services for obtaining characteristic matrices
 * for several material modes is provided depending on other abstract services.
 * Also general implementation of giveRealStressVector service is provided,
 * computing the stress increment vector from strain increment multiplied by
 * stiffness.
 *
 * Tasks:
 * - Returning standard material stiffness and flexibility matrices for 3d-case.
 *   according to current state determined by using data stored in Gausspoint.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at gauss point for 3d - case.
 */
class LinearElasticMaterial : public StructuralMaterial
{
public:
    /// Constructor.
    LinearElasticMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }
    /// Destructor.
    virtual ~LinearElasticMaterial() { }

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strain, TimeStep *tStep);
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStressVector_3dDegeneratedShell(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep);
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStressVector_Warping(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);
    virtual void giveRealStressVector_2dBeamLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    virtual void giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    virtual void giveRealStressVector_Fiber(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);

    virtual void giveEshelbyStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep);
    double giveEnergyDensity(GaussPoint *gp, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual double giveShearModulus() { return 1.; }
    virtual int hasNonLinearBehaviour() { return 0; }
    virtual const char *giveClassName() const { return "LinearElasticMaterial"; }
};
} // end namespace oofem
#endif // linearelasticmaterial_h
