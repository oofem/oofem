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

#ifndef levelsetpcs_h
#define levelsetpcs_h

#include "materialinterface.h"
#include "mathfem.h"
#include "geotoolbox.h"
#ifndef __MAKEDEPEND
 #include <vector>
#endif

namespace oofem {
#define LevelSetPCS_CACHE_ELEMENT_VOF 0

class LevelSetPCS;

/**
 * Element interface for LevelSetPCS class representing level-set like material interface.
 * The elements should provide specific functionality in order to colaborate with LEPlic and this
 * required functionality is declared in this interface.
 */
class LevelSetPCSElementInterface : public Interface
{
public:
    LevelSetPCSElementInterface() { }

    /// @name The element interface required by LevelSetPCSElementInterface.
    //@{
    /**
     * Evaluates F in level set equation of the form
     * @f[ \phi_t + F(\nabla\phi, x) |\nabla\phi| = 0 @f]
     * where for interface position driven by flow with speed u:
     * @f[ F = u\cdot \frac{\nabla\phi}{|\nabla\phi|} @f]
     */
    virtual double LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep) = 0;

    /**
     * Returns gradient of shape functions.
     */
    virtual void LS_PCS_computedN(FloatMatrix &answer) = 0;
    /// Returns receiver's volume
    virtual double LS_PCS_computeVolume() = 0;

    /**
     * Evaluates S in level set equation of the form
     * @f[ \phi_t = S(\phi) (1-|\nabla\phi|) = 0 @f]
     * where
     * @f[ S=\frac{\phi}{\sqrt{\phi^2+\epsilon^2}} @f]
     */
    virtual double LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep) = 0;

    /**
     * Returns VOF fractions for each material on element
     * according to nodal values of level set function (passed as parameter)
     */
    virtual void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi) = 0;
    //@}
};

/**
 * Abstract base class representing Level Set representation of material interfaces.
 * The solution algorithm is based on positive coefficint scheme.
 * This algorithm is limitted to d-dimensional simplexes with linear approximation.
 * Its typical use to model moving interface (such as free surface)
 * in a fixed-grid methods (as typically used in CFD).
 * The basic tasks are representation of interface and its updating.
 */
class LevelSetPCS : public MaterialInterface
{
protected:
    /// Array used to store value of level set function for each node.
    FloatArray levelSetValues, previousLevelSetValues;
    enum PCSEqType { PCS_levelSetUpdate, PCS_levelSetRedistance };
    Polygon initialRefMatVol;
    bool initialRefMatFlag;
    /// Indexes of nodal coordinates used to init levelset using initialRefMatVol.
    int ci1, ci2;

    /// Type of reinitialization algorithm to use.
    int reinit_alg;
    /// Time step used in reinitialization of LS (if apply).
    double reinit_dt;
    bool reinit_dt_flag;
    /// Reinitialization error limit.
    double reinit_err;
    /// number of spatial dimensions.
    int nsd;
    /// Level set values version.
    long int levelSetVersion;

#ifdef LevelSetPCS_CACHE_ELEMENT_VOF
    std :: vector< FloatArray >elemVof;
    long int elemVofLevelSetVersion;
#endif

public:
    /** Constructor. Takes two two arguments. Creates
     *  MaterialInterface instance with given number and belonging to given domain.
     *  @param n Component number in particular domain. For instance, can represent
     *  node number in particular domain.
     *  @param d Domain to which component belongs to.
     */
    LevelSetPCS(int n, Domain *d) : MaterialInterface(n, d) {
        initialRefMatFlag = false;
        reinit_dt_flag = false;
        levelSetVersion = 0;
    }

    virtual void initialize();
    virtual void updatePosition(TimeStep *atTime);
    virtual void updateYourself(TimeStep *tStep) { previousLevelSetValues = levelSetValues; }
    double computeCriticalTimeStep(TimeStep *);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void reinitialization(TimeStep *atTime);

    virtual void giveMaterialMixtureAt(FloatArray &answer, FloatArray &position);
    virtual void giveElementMaterialMixture(FloatArray &answer, int ielem);
    virtual double giveNodalScalarRepresentation(int i) { return levelSetValues.at(i); }

    /// Returns level set value in specific node
    double giveLevelSetDofManValue(int i) { return levelSetValues.at(i); }

    // identification
    const char *giveClassName() const { return "LevelSetPCS"; }
    classType giveClassID() const { return LevelSetPCSClass; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

protected:
    void pcs_stage1(FloatArray &ls, FloatArray &fs, FloatArray &w, TimeStep *atTime, PCSEqType t);
    double evalElemFContribution(PCSEqType t, int ie, TimeStep *atTime);
    double evalElemfContribution(PCSEqType t, int ie, TimeStep *atTime);

    /**
     * Reinitializes the level set representation by solving
     * @f$ d_{\tau} = S(\phi)(1-|\nabla d|) @f$ to steady state.
     */
    void redistance(TimeStep *atTime);


    /** @name Fast marching related services */
    //@{
    /** Reinitializes the level set representation using fast marching method. */
    void FMMReinitialization(FloatArray &ls);
    //@}
};
} // end namespace oofem
#endif // levelsetpcs_h
