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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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
 * Class representing material status for the mesoscale sintering material
 * @author Mikael Öhman
 */
class FE2SinteringMaterialStatus : public StructuralMaterialStatus
{
protected:
    double volume;

    // TopologyDescription td;

    /// Input file (for now)
    char file[1024];

    StokesFlowStressHomogenization *rve;

	FloatArray oldStrainVector;

public:

    /** Creates new material status.
     * @param n number ?
     * @param d domain that status belongs to
     * @param g gauss point that the status belongs to
     */
	FE2SinteringMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~FE2SinteringMaterialStatus() { delete this->rve; }

    double giveVolume() { return this->volume; }
    StokesFlowStressHomogenization *giveRVE() { return this->rve; }
    //HomogenizationInterface *giveHomogenizationInterface() { return this->hrve; }

    FloatArray &giveOldStrainVector() { return this->oldStrainVector; }

    /// Creates/Iniciates the RVE problem
    bool createRVE(int n, GaussPoint *gp);

    /// Copies time step data to RVE
    void setTimeStep(TimeStep *tStep);

    /// Print receiver's output to given stream.
    void printOutputAt(FILE *, TimeStep *);

    /**
     * Initializes the temporary internal variables (stresss and strains vectors),
     * describing the current state according to
     * previously reached equilibrium internal variables
     */
    virtual void initTempStatus();
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *);

    /**
     * Stores context of receiver into given stream (the equilibriun stress and strains vectors are stored).
     * Generally, only non-temp internal history variables should be stored.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream (the equilibriun stress and strains vectors are restored).
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType
     * @exception throws an ContextIOERR exception if error encountered
     */
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
protected:
    /// viscosity
    double mu;
    /// surface energy
    double gamma_s;

public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    FE2SinteringMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }
    /// Destructor.
    ~FE2SinteringMaterial() { }

    virtual void giveRealStressVector (FloatArray &answer, MatResponseForm, GaussPoint *, const FloatArray &, TimeStep *);
    virtual void givePlaneStressStiffMtrx (FloatMatrix &answer, MatResponseForm, MatResponseMode, GaussPoint *gp, TimeStep *atTime);

    /**
     * Initializes receiver acording to object description stored in input record.
     * The density of material is read into property dictionary (keyword 'd')
     */
    IRResultType initializeFrom(InputRecord *ir);

    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    /**
     * Tests, if material supports material mode.
     * @param mode required material mode
     * @return nonzero if supported, zero otherwise
     */
    virtual int hasMaterialModeCapability(MaterialMode mode);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "FE2SinteringMaterial"; }

    /// Returns classType id of receiver.
    classType giveClassID() const { return FE2SinteringMaterialClass; }

    /** Allows programmer to test some internal data, before computation begins.
     *  For example, one may use this function, to ensure that element has material with
     *  required capabilities is assigned to element. This must be done after all
     *  mesh components are instanciated.
     *  @return nonzero if receiver check is o.k.
     */
    virtual int checkConsistency();

    /**
     * Creates new copy of associated status and inserts it into given integration point.
     * @param gp Integration point where newly created status will be stored.
     * @return reference to new status.
     */
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

private:
    static int n;
};

} // end namespace oofem
#endif // rvesinteringmaterial_h
