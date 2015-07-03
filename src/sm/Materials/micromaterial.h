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

#ifndef micromaterial_h
#define micromaterial_h

#include "structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "datastream.h"
#include "contextioerr.h"
#include "unknownnumberingscheme.h"
#include "boundarycondition.h"
#include "Elements/3D/macrolspace.h"
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
    MicroMaterialStatus(int, Domain * d, GaussPoint * gp);

    /// Destructor
    virtual ~MicroMaterialStatus();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "MicroMaterialStatus"; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
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
    /// Destructor
    virtual ~MicroMaterial();

    std :: string inputFileNameMicro;

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveInputRecordName() const { return _IFT_MicroMaterial_Name; }
    virtual const char *giveClassName() const { return "MicroMaterial"; }

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *, const FloatArray &, TimeStep *);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    void giveMacroStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep, MatResponseMode rMode, const IntArray &microMasterNodes, const IntArray &microBoundaryNodes);

    void setMacroProperties(Domain *macroDomain, MacroLSpace *macroLSpaceElement, const IntArray &microMasterNodes, const IntArray &microBoundaryNodes);

    /// Pointer to the underlying micro problem.
    EngngModel *problemMicro;

    /// Pointer to the macroscale domain.
    Domain *macroDomain;

    /// Pointer to the macroscale element.
    MacroLSpace *macroLSpaceElement;

    /// Related to numbering scheme.
    void init(void);
    int giveDofEquationNumber(Dof *dof) const;
    virtual bool isDefault() const { return isDefaultNumbering; }
    virtual int giveRequiredNumberOfDomainEquation() const;
    //friend class EngngModel;-not here but define in EngngModel class
    /// Array containing coordinates of 8 master nodes of microproblem.
    std::vector< FloatArray >microMasterCoords;
    /// Array containing equation numbers for boundary nodes [DofManagerNumber][DOF].
    int **microBoundaryDofs;
    /// Array of equation numbers associated to boundary nodes.
    IntArray microBoundaryDofsArr;
    /// Array containing equation numbers for internal nodes to be condensed out [DofManagerNumber][DOF].
    int **microInternalDofs;
    /// Array of equation numbers associated to internal nodes.
    IntArray microInternalDofsArr;
    /// Array containing default equation numbers for all nodes [DofManagerNumber][DOF].
    int **microDefaultDofs;
    /// Flag signalizing whether micromaterial is used by other element.
    bool microMatIsUsed;

protected:
    bool isDefaultNumbering;
    /// The maximum DOFs corresponding to released all of the boundary conditions.
    int maxNumberOfDomainEquation;
    /// Required number of domain equations.
    int reqNumberOfDomainEquation;
    /// Number of DOF Managers.
    int NumberOfDofManagers;
    enum EquationNumbering { AllNodes, BoundaryNodes, InteriorNodes };
    EquationNumbering DofEquationNumbering;
    /// Number of equations associated with boundary nodes.
    int totalBoundaryDofs;
    /// Number of equations associated with boundary nodes.
    int totalInternalDofs;
};
} // end namespace oofem
#endif // micromaterial_h
