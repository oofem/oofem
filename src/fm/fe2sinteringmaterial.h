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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef fe2sinteringmaterial_h
#define fe2sinteringmaterial_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "structuralelement.h"
#include "matstatus.h"
#include "stokesflowstresshomogenization.h"

namespace oofem {

class GaussPoint;

/**
 * Class representing material status for the mesoscale sintering material.
 * @author Mikael Öhman
 */
class FE2SinteringMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Input file
    char file[OOFEM_MAX_LINE_LENGTH];

    StokesFlowStressHomogenization *rve;

	double voffraction;

public:
    /**
     * Creates new material status.
     * @param n Material status number.
     * @param d Domain that status belongs to.
     * @param gp Gauss point that the status belongs to.
     * @param inputfile The input file describing the micro problem.
     * @param porosity Initial porosity.
     */
	FE2SinteringMaterialStatus(int n, Domain *d, GaussPoint *gp, const char* inputfile);
    /// Destructor
    ~FE2SinteringMaterialStatus() { delete this->rve; }

    StokesFlowStressHomogenization *giveRVE() { return this->rve; }

    double giveVOFFraction() { return this->voffraction; }

    /// Creates/Initiates the RVE problem.
    bool createRVE(int n, GaussPoint *gp, const char *inputfile);

    /// Copies time step data to RVE.
    void setTimeStep(TimeStep *tStep);

    void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    const char *giveClassName() const { return "FE2SinteringMaterialStatus"; }
    classType giveClassID() const { return FE2SinteringMaterialStatusClass; }
};


/**
 * Multiscale constitutive model of free surface tension driven fluids (at mesoscale)
 * @author Mikael Öhman
 */
class FE2SinteringMaterial : public StructuralMaterial
{
private:
    char inputfile[OOFEM_MAX_LINE_LENGTH];

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    FE2SinteringMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }
    /// Destructor.
    ~FE2SinteringMaterial() { }

    virtual void giveRealStressVector (FloatArray &answer, MatResponseForm, GaussPoint *, const FloatArray &, TimeStep *);
    virtual void givePlaneStressStiffMtrx (FloatMatrix &answer, MatResponseForm, MatResponseMode, GaussPoint *gp, TimeStep *atTime);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "FE2SinteringMaterial"; }
    classType giveClassID() const { return FE2SinteringMaterialClass; }
    virtual int checkConsistency();
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

private:
    static int n;
};

} // end namespace oofem
#endif // rvesinteringmaterial_h
