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

#ifndef SRC_SM_MATERIALS_FE2STRUCTURALMATERIAL_H_
#define SRC_SM_MATERIALS_FE2STRUCTURALMATERIAL_H_

#include "structuralmaterial.h"
#include "prescribedgradientbcweak.h"
#include "structuralms.h"

#include <memory>


///@name Input fields for FE^2 structural material
//@{
#define _IFT_FE2StructuralMaterial_Name "fe2structuralmaterial"
#define _IFT_FE2StructuralMaterial_fileName "inputfile"
//@}

namespace oofem {
class EngngModel;

/**
 * Class representing material status for the subscale structure, i.e a Representative Volume Element (RVE).
 * @author Erik Svenning
 */
class FE2StructuralMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// The subscale microstructure
    std :: unique_ptr< EngngModel > rve;

    /// Boundary condition on RVE that performs the computational homogenization.
    PrescribedGradientBCWeak *bc;

public:
    /**
     * Creates new material status.
     * @param n Material status number.
     * @param d Domain that status belongs to.
     * @param gp Gauss point that the status belongs to.
     * @param inputfile The input file describing the micro problem.
     */
    FE2StructuralMaterialStatus(int n, Domain * d, GaussPoint * gp, const std :: string & inputfile);
    /// Destructor
    virtual ~FE2StructuralMaterialStatus();

    EngngModel *giveRVE() { return this->rve.get(); }
    PrescribedGradientBCWeak *giveBC() { return this->bc; }

    /// Initiates the RVE problem.
    bool createRVE(int n, GaussPoint *gp, const std :: string &inputfile);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "FE2StructuralMaterialStatus"; }
};

/**
 * Multiscale constitutive model representing the microstructure of the material as a Representative Volume Element (RVE).
 * @author Erik Svenning
 */
class FE2StructuralMaterial : public StructuralMaterial {
public:
	FE2StructuralMaterial(int n, Domain * d);
	virtual ~FE2StructuralMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual int checkConsistency() {return true;}

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual const char *giveClassName() const { return "FE2StructuralMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_FE2StructuralMaterial_Name; }

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);

    /**
     * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
     * in integration point.
     * @param answer Computed results.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

private:
    std :: string inputfile;
    static int n;
};

} /* namespace oofem */

#endif /* SRC_SM_MATERIALS_FE2STRUCTURALMATERIAL_H_ */
