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

#ifndef levelsetpcs_h
#define levelsetpcs_h

#include "materialinterface.h"
#include "geotoolbox.h"
#include "interface.h"

#include <vector>

///@name Input fields for LSPCS
//@{
#define _IFT_LevelSetPCS_levelSetValues "levelset"
#define _IFT_LevelSetPCS_refmatpoly_x "refmatpolyx"
#define _IFT_LevelSetPCS_refmatpoly_y "refmatpolyy"
#define _IFT_LevelSetPCS_reinit_dt "rdt"
#define _IFT_LevelSetPCS_reinit_err "rerr"
#define _IFT_LevelSetPCS_reinit_alg "lsra"
#define _IFT_LevelSetPCS_nsd "nsd"
#define _IFT_LevelSetPCS_ci1 "ci1"
#define _IFT_LevelSetPCS_ci2 "ci2"
//@}


namespace oofem {
#define LevelSetPCS_CACHE_ELEMENT_VOF 0

class LevelSetPCS;

/**
 * Element interface for LevelSetPCS class representing level-set like material interface.
 * The elements should provide specific functionality in order to collaborate with LEPlic and this
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
 * The solution algorithm is based on positive coefficient scheme.
 * This algorithm is limited to d-dimensional simplexes with linear approximation.
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
     *  a MaterialInterface instance with given number and belonging to given domain.
     *  @param n Component number in particular domain. For instance, can represent
     *  node number in particular domain.
     *  @param d Domain to which component belongs to.
     */
    LevelSetPCS(int n, Domain * d) : MaterialInterface(n, d) {
        initialRefMatFlag = false;
        reinit_dt_flag = false;
        levelSetVersion = 0;
    }

    virtual void initialize();
    virtual void updatePosition(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep) { previousLevelSetValues = levelSetValues; }
    virtual double computeCriticalTimeStep(TimeStep *tStep);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual void reinitialization(TimeStep *tStep);

    virtual void giveMaterialMixtureAt(FloatArray &answer, FloatArray &position);
    virtual void giveElementMaterialMixture(FloatArray &answer, int ielem);
    virtual double giveNodalScalarRepresentation(int i) { return levelSetValues.at(i); }

    /// Returns level set value in specific node
    double giveLevelSetDofManValue(int i) { return levelSetValues.at(i); }

    // identification
    virtual const char *giveClassName() const { return "LevelSetPCS"; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

protected:
    void pcs_stage1(FloatArray &ls, FloatArray &fs, FloatArray &w, TimeStep *tStep, PCSEqType t);
    double evalElemFContribution(PCSEqType t, int ie, TimeStep *tStep);
    double evalElemfContribution(PCSEqType t, int ie, TimeStep *tStep);

    /**
     * Reinitializes the level set representation by solving
     * @f$ d_{\tau} = S(\phi)(1-|\nabla d|) @f$ to steady state.
     */
    void redistance(TimeStep *tStep);

    /** @name Fast marching related services */
    //@{
    /** Reinitializes the level set representation using fast marching method. */
    void FMMReinitialization(FloatArray &ls);
    //@}
};
} // end namespace oofem
#endif // levelsetpcs_h
