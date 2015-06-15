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

#ifndef pfem_h
#define pfem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "dofdistributedprimaryfield.h"
//<RESTRICTED_SECTION>
#include "materialinterface.h"
//</RESTRICTED_SECTION>
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#include "pfemnumberingschemes.h"

///@name Input fields for PFEM
//@{
#define _IFT_PFEM_Name "pfem"
#define _IFT_PFEM_deltat "deltat"
#define _IFT_PFEM_mindeltat "mindeltat"
#define _IFT_PFEM_alphashapecoef "alphashapecoef"
#define _IFT_PFEM_particalRemovalRatio "removalratio"
#define _IFT_PFEM_printVolumeReport "volumereport"
#define _IFT_PFEM_discretizationScheme "scheme"
#define _IFT_PFEM_associatedMaterial "material"
#define _IFT_PFEM_associatedCrossSection "cs"
#define _IFT_PFEM_pressureBC "pressure"

//@}

namespace oofem {

class PressureRhsAssembler : public VectorAssembler
{
public:
    virtual void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const;
};


/**
 * This class represents PFEM method for solving incompressible Navier-Stokes equations
 */
class PFEM : public EngngModel
{
protected:
    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

    /// left-hand side matrix for the auxiliary velocity equations
    FloatArray avLhs;
    /// left-hand side matrix for the pressure equations
    SparseMtrx *pLhs;
    /// left-hand side matrix for the velocity equations
    FloatArray vLhs;

    // TODO: consider using DofDistributedPrimaryField for pressure and velocity
    /// Pressure field
    PrimaryField PressureField;
    /// Velocity field
    PrimaryField VelocityField;

    FloatArray AuxVelocity;
    FloatArray prescribedTractionPressure;

    /// time step
    double deltaT;
    /// minimal value of time step
    double minDeltaT;
    /// value of alpha coefficient for the boundary recognition
    double alphaShapeCoef;
    /// element side ratio for the removal of the close partices
    double particleRemovalRatio;
    
    /// area or volume of the fluid domain, which can be controlled
    double domainVolume;
    /// flag for volume report
    bool printVolumeReport;

    /// explicit or implicit time discretization
    int discretizationScheme;

    /// number of cross section to associate with created elements
    int associatedCrossSection;
    /// number of material to associate with created elements
    int associatedMaterial;
    /// number of pressure boundary condition to be prescribed on free surface
    int associatedPressureBC;

    /// Pressure numbering
    PressureNumberingScheme pns;
    /// Auxiliary Velocity numbering
    AuxVelocityNumberingScheme avns;
    /// Velocity numbering
    VelocityNumberingScheme vns;
    /// Prescribed velocity numbering
    VelocityNumberingScheme prescribedVns;

public:
    PFEM(int i, EngngModel * _master = NULL) :
        EngngModel(i, _master)
        , avLhs()
        , pLhs(NULL)
        , vLhs()
        , PressureField(this, 1, FT_Pressure, 1)
        , VelocityField(this, 1, FT_Velocity, 1)
        , pns()
        , vns(false)
        , prescribedVns(true)
    {
        ndomains = 1;
        nMethod = NULL;
        domainVolume = 0.0;
        printVolumeReport = false;
        discretizationScheme = 1; // implicit iterative scheme is default
        associatedCrossSection = 0;
        associatedMaterial = 0;
        associatedPressureBC = 0;
    }
    ~PFEM() { }

    void solveYourselfAt(TimeStep *);
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.giv
     */
    virtual void               updateYourself(TimeStep *);

    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);


    TimeStep *giveNextStep();
    TimeStep *giveSolutionStepWhenIcApply();
    NumericalMethod *giveNumericalMethod(MetaStep *);

    /** Removes all elements and call DelaunayTriangulator to build up new mesh with new recognized boundary.
     *  Non-meshed particles are set free and move according Newton laws of motion
     */
    virtual void preInitializeNextStep();

    //equation numbering using PressureNumberingScheme
    virtual int forceEquationNumbering(int id);
    
    virtual int requiresUnknownsDictionaryUpdate() { return true; }


    /// Initialization from given input record
    IRResultType initializeFrom(InputRecord *ir);

    // consistency check
    virtual int checkConsistency(); // returns nonzero if o.k.
    // identification
    const char *giveClassName() const { return "PFEM"; }
    fMode giveFormulation() { return AL; }
    //    fMode giveFormulation() { return TL; }

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     *  Dof class provides component printing routines, but emodel is responsible
     *  for what will be printed at DOF level.
     *  @param stream output stream
     *  @param iDof dof to be processed
     *  @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    virtual int giveNumberOfDomainEquations(int, const UnknownNumberingScheme &num);

    virtual int giveNewEquationNumber(int domain, DofIDItem);
    virtual int giveNewPrescribedEquationNumber(int domain, DofIDItem);

    virtual int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *stepN); // { return ( int ) mode; }
    virtual void updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep);

    /// Writes pressures into the dof unknown dictionaries
    void updateDofUnknownsDictionaryPressure(DofManager *inode, TimeStep *tStep);
    /// Writes velocities into the dof unknown dictionaries
    void updateDofUnknownsDictionaryVelocities(DofManager *inode, TimeStep *tStep);
    /// Resets the equation numberings as the mesh is recreated in each time step
    void resetEquationNumberings();
    
    /// Returns number of material to be associated with created elements
    int giveAssociatedMaterialNumber() {return associatedMaterial;}
    /// Returns number of cross section to be associated with created elements
    int giveAssociatedCrossSectionNumber() {return associatedCrossSection;}
    /// Returns number of zero pressure boundary condition to be prescribed on free surface
    int giveAssociatedPressureBC() {return associatedPressureBC;}

protected:
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     */
    void updateInternalState(TimeStep *);
    /// Initializes velocity and pressure fields in the first step with prescribed values
    void applyIC(TimeStep *);
    /// Deactivates particles upon the particalRemovalRatio
    void deactivateTooCloseParticles();
};
} // end namespace oofem
#endif // pfem_h
