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

#ifndef structuralfe2material_h
#define structuralfe2material_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"

#include <memory>

///@name Input fields for StructuralFE2Material
//@{
#define _IFT_StructuralFE2Material_Name "structfe2material"
#define _IFT_StructuralFE2Material_fileName "filename"
#define _IFT_StructuralFE2Material_useNumericalTangent "use_num_tangent"
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
    PrescribedGradientHomogenization *bc = nullptr;

    FloatMatrix tangent;
    bool oldTangent = true;

    /// Interface normal direction
    FloatArray mNormalDir;

    std :: string mInputFile;

public:
    StructuralFE2MaterialStatus(int rank, GaussPoint * g,  const std :: string & inputfile);

    EngngModel *giveRVE() const { return this->rve.get(); }
    PrescribedGradientHomogenization *giveBC();// { return this->bc; }

    void markOldTangent();
    void computeTangent(TimeStep *tStep);

    /// Creates/Initiates the RVE problem.
    bool createRVE(int n, const std :: string &inputfile, int rank);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    FloatMatrix &giveTangent() { return tangent; }

    const char *giveClassName() const override { return "StructuralFE2MaterialStatus"; }

    void initTempStatus() override;

    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const FloatArray &giveNormal() const { return mNormalDir; }
    void letNormalBe(FloatArray iN) { mNormalDir = std :: move(iN); }

    double giveRveLength();

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables(const MaterialStatus &iStatus) override;
    void addStateVariables(const MaterialStatus &iStatus) override { OOFEM_ERROR("Not implemented."); }

    // For debugging only
    bool mNewlyInitialized = true;
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
    bool useNumTangent = false;

public:
    StructuralFE2Material(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveInputRecordName() const override { return _IFT_StructuralFE2Material_Name; }
    const char *giveClassName() const override { return "StructuralFE2Material"; }
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return true; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
};

} // end namespace oofem
#endif // structuralfe2material_h
