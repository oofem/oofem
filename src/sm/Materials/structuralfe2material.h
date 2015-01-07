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

#ifndef tutorialmaterial_h
#define tutorialmaterial_h

#include "Materials/structuralmaterial.h"
#include "Materials/structuralms.h"

#include <memory>

///@name Input fields for StructuralFE2Material
//@{
#define _IFT_StructuralFE2Material_Name "structfe2material"
#define _IFT_StructuralFE2Material_fileName "filename"
//@}

namespace oofem {
class EngngModel;
class PrescribedGradientHomogenization;

class StructuralFE2MaterialStatus : public StructuralMaterialStatus
{
protected:
    /// The RVE
    std :: unique_ptr< EngngModel > rve;
    /// Boundary condition in RVE that performs the computational homogenization.
    PrescribedGradientHomogenization *bc;

    FloatMatrix tangent;
    bool oldTangent;

public:
    StructuralFE2MaterialStatus(int n, Domain * d, GaussPoint * g,  const std :: string & inputfile);
    virtual ~StructuralFE2MaterialStatus() {}

    EngngModel *giveRVE() { return this->rve.get(); }
    PrescribedGradientHomogenization *giveBC() { return this->bc; }

    void markOldTangent();
    void computeTangent(TimeStep *tStep);

    /// Creates/Initiates the RVE problem.
    bool createRVE(int n, GaussPoint *gp, const std :: string &inputfile);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    FloatMatrix &giveTangent() { return tangent; }
    
    virtual const char *giveClassName() const { return "StructuralFE2MaterialStatus"; }
    
    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};

    
/**
 * Multiscale constitutive model for subscale structural problems.
 *
 * The material uses the PrescribedGradient boundary condition to perform computational homogenization.
 * The requirement for the supplied subscale problem is:
 * - It must have a PrescribedGradient boundary condition.
 * - It must be the first boundary condition
 *
 * @author Mikael Öhman 
 */
class StructuralFE2Material : public StructuralMaterial
{
protected:
    std :: string inputfile;
    static int n;

public:
    StructuralFE2Material(int n, Domain * d);
    virtual ~StructuralFE2Material();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveInputRecordName() const { return _IFT_StructuralFE2Material_Name; }
    virtual const char *giveClassName() const { return "StructuralFE2Material"; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    const void giveDeviatoricProjectionMatrix(FloatMatrix &answer);
    // stress computation methods
    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
};

} // end namespace oofem
#endif // structuralfe2material
