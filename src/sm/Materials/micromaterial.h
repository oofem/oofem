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

#ifndef micromaterial_h
#define micromaterial_h

#include "structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "datastream.h"
#include "contextioerr.h"
#include "unknownnumberingscheme.h"
#include "boundarycondition.h"
#include "sm/Elements/3D/macrolspace.h"
#include "error.h"

///@name Input fields for MicroMaterial
//@{
#define _IFT_MicroMaterial_Name "micromat"
#define _IFT_MicroMaterial_fileName "file"
//@}

namespace oofem {
class UnknownNumberingScheme;
class MicroMaterial;
class MacroLSpace;

class MicroMaterialStatus : public StructuralMaterialStatus
{
public:
    /// Constructor
    MicroMaterialStatus(GaussPoint * gp);

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "MicroMaterialStatus"; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * This class is a base class for microproblem.
 * The microproblem represents itself a problem which is solved separately from the macroproblem with appropriate boundary conditions.
 * MacroLspace needs stiffness matrix derived from this microproblem.
 * For this purpose, natural boundary conditions on microproblem have to be excluded.
 * Stiffness matrix of microproblem is condensed to provide stiffness matrix for macroelement.
 */
class MicroMaterial : public StructuralMaterial, public UnknownNumberingScheme
{
public:
    /// Constructor
    MicroMaterial(int n, Domain * d);

    std :: string inputFileNameMicro;

    void initializeFrom(InputRecord &ir) override;

    const char *giveInputRecordName() const override { return _IFT_MicroMaterial_Name; }
    const char *giveClassName() const override { return "MicroMaterial"; }

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    void giveMacroStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep, MatResponseMode rMode, const IntArray &microMasterNodes, const IntArray &microBoundaryNodes);

    void setMacroProperties(Domain *macroDomain, MacroLSpace *macroLSpaceElement, const IntArray &microMasterNodes, const IntArray &microBoundaryNodes);

    /// Pointer to the underlying micro problem.
    std::unique_ptr<EngngModel> problemMicro;

    /// Pointer to the macroscale domain.
    Domain *macroDomain = nullptr;

    /// Pointer to the macroscale element.
    MacroLSpace *macroLSpaceElement = nullptr;

    /// Related to numbering scheme.
    void init(void) override;
    int giveDofEquationNumber(Dof *dof) const override;
    bool isDefault() const override { return isDefaultNumbering; }
    int giveRequiredNumberOfDomainEquation() const override;
    //friend class EngngModel;-not here but define in EngngModel class
    /// Array containing coordinates of 8 master nodes of microproblem.
    std::vector< FloatArray >microMasterCoords;
    /// Array containing equation numbers for boundary nodes [DofManagerNumber][DOF].
    std::vector<IntArray> microBoundaryDofs;
    /// Array of equation numbers associated to boundary nodes.
    IntArray microBoundaryDofsArr;
    /// Array containing equation numbers for internal nodes to be condensed out [DofManagerNumber][DOF].
    std::vector<IntArray> microInternalDofs;
    /// Array of equation numbers associated to internal nodes.
    IntArray microInternalDofsArr;
    /// Array containing default equation numbers for all nodes [DofManagerNumber][DOF].
    std::vector<IntArray> microDefaultDofs;
    /// Flag signalizing whether micromaterial is used by other element.
    bool microMatIsUsed = false;

protected:
    bool isDefaultNumbering = true;
    /// The maximum DOFs corresponding to released all of the boundary conditions.
    int maxNumberOfDomainEquation = 0;
    /// Required number of domain equations.
    int reqNumberOfDomainEquation = 0;
    /// Number of DOF Managers.
    int NumberOfDofManagers = 0;
    enum EquationNumbering { AllNodes, BoundaryNodes, InteriorNodes };
    EquationNumbering DofEquationNumbering = AllNodes;
    /// Number of equations associated with boundary nodes.
    int totalBoundaryDofs = 0;
    /// Number of equations associated with boundary nodes.
    int totalInternalDofs = 0;
};
} // end namespace oofem
#endif // micromaterial_h
