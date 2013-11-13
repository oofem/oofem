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

#ifndef tr_shell01_h
#define tr_shell01_h

#include "structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "zzerrorestimator.h"
#include "cct3d.h"
#include "trplanrot3d.h"
#include "spatiallocalizer.h"

#define _IFT_TR_SHELL01_Name "tr_shell01"

namespace oofem {

/**
 * This class implements an triangular three-node shell finite element, composed of
 * cct3d and trplanrot3d elements.
 * Each node has 6 degrees of freedom.
 *
 * @author Ladislav Svoboda
 */
class TR_SHELL01 : public StructuralElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface, public SpatialLocalizerInterface
{
protected:
    /// Pointer to plate element.
    CCTPlate3d *plate;
    /// Pointer to membrane (plane stress) element.
    TrPlaneStrRot3d *membrane;
    /**
     * Element integraton rule (plate and membrane parts have their own integration rules)
     * this one used to integrate element error and perhaps can be (re)used for other putrposes.
     * Created on demand.
     */
    IntegrationRule *compositeIR;

public:
    /// Constructor
    TR_SHELL01(int n, Domain *d);
    /// Destructor
    virtual ~TR_SHELL01() {
        delete plate;
        delete membrane;
        if (this->compositeIR) delete this->compositeIR;
    }

    virtual FEInterpolation *giveInterpolation() const { return plate->giveInterpolation(); }

    virtual int computeNumberOfDofs() { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
    { plate->giveDofManDofIDMask(inode, ut, answer); }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_TR_SHELL01_Name; }
    virtual const char *giveClassName() const { return "TR_SHELL01"; }
    virtual classType giveClassID() const { return TR_SHELL01Class; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual bool giveRotationMatrix(FloatMatrix &answer, EquationID eid);

    virtual void updateYourself(TimeStep *tStep);
    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

#ifdef __OOFEG
    virtual void drawRawGeometry(oofegGraphicContext &);
    virtual void drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &context);
#endif
    // the membrane and plate irules are same (chacked in initializeFrom)
    virtual int giveDefaultIntegrationRule() const { return plate->giveDefaultIntegrationRule();}
    virtual IntegrationRule *giveDefaultIntegrationRulePtr() {return plate->giveDefaultIntegrationRulePtr();}
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
    virtual integrationDomain giveIntegrationDomain() const { return _Triangle; }
    virtual MaterialMode giveMaterialMode() { return _Unknown; }

    virtual Interface *giveInterface(InterfaceType it);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    // ZZErrorEstimatorInterface
    virtual Element *ZZErrorEstimatorI_giveElement() { return this; }

    virtual IntegrationRule *ZZErrorEstimatorI_giveIntegrationRule();
    virtual void ZZErrorEstimatorI_computeLocalStress(FloatArray& answer, FloatArray& sig);

    // ZZRemeshingCriteriaInterface
    virtual double ZZRemeshingCriteriaI_giveCharacteristicSize();
    virtual int ZZRemeshingCriteriaI_givePolynOrder() { return 1; };

    // SpatialLocalizerI
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);
    virtual void SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1);


    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) {
        return this->plate->computeGlobalCoordinates (answer, lcoords);
    }

protected:
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS)
    { _error("TR_SHELL01 :: computeBmatrixAt: calling of this function is not allowed"); }
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &)
    { _error("TR_SHELL01 :: computeNmatrixAt: calling of this function is not allowed"); }

    /// @todo In time delete
protected:
    virtual void computeGaussPoints()
    { this->membrane->computeGaussPoints(); this->plate->computeGaussPoints(); }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *stepN)
    { _error("TR_SHELL01 :: computeStressVector: calling of this function is not allowed"); }
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }

public:
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
    { _error("TR_SHELL01 :: ...: calling of this function is not allowed"); }
};
} // end namespace oofem
#endif
