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

#ifndef tr_shell01_h
#define tr_shell01_h

#include "structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "zzerrorestimator.h"
#include "cct3d.h"
#include "trplanrot3d.h"

namespace oofem {

/**
 * This class implements an triangular three-node shell finite element, composed of
 * cct3d and trplanrot3d elements.
 * Each node has 6 degrees of freedom.
 *
 * @author Ladislav Svoboda
 */
  class TR_SHELL01 : public StructuralElement, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, 
    public ZZErrorEstimatorInterface, public ZZRemeshingCriteriaInterface
{
protected:
    /// Pointer to plate element.
    CCTPlate3d *plate;
    /// Pointer to membrane (plane stress) element.
    TrPlaneStrRot3d *membrane;

public:
    /// Constructor
    TR_SHELL01(int n, Domain *d);
    /// Destructor
    virtual ~TR_SHELL01() {
        delete plate;
        delete membrane;
    }

    virtual int computeNumberOfDofs(EquationID ut) { return 18; }
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
    { plate->giveDofManDofIDMask(inode, ut, answer); }
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);
    // definition & identification
    virtual const char *giveClassName() const { return "TR_SHELL01"; }
    virtual classType giveClassID() const { return TR_SHELL01Class; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep);
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void updateYourself(TimeStep *tStep);
    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType type);
    virtual void drawScalar(oofegGraphicContext &context);
    //void drawInternalState(oofegGraphicContext &);
#endif
    // the membrane and plate irules are same (chacked in initializeFrom)
    virtual int giveDefaultIntegrationRule() const { return plate->giveDefaultIntegrationRule();}
    virtual IntegrationRule *giveDefaultIntegrationRulePtr() {return plate->giveDefaultIntegrationRulePtr();}
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
    virtual integrationDomain giveIntegrationDomain() { return _Triangle; }
    virtual MaterialMode giveMaterialMode() { return _Unknown; }

    virtual Interface *giveInterface(InterfaceType it);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }
    virtual void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type);

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    // ZZErrorEstimatorInterface
    virtual Element *ZZErrorEstimatorI_giveElement() { return this; }
    virtual void ZZErrorEstimatorI_computeEstimatedStressInterpolationMtrx(FloatArray &answer, GaussPoint *gp,
                                                                           InternalStateType type)
    { ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(answer, gp, type); }

    // ZZRemeshingCriteriaInterface
    virtual double ZZRemeshingCriteriaI_giveCharacteristicSize();
    virtual int ZZRemeshingCriteriaI_givePolynOrder() { return 1; };


protected:
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS)
    { _error("TR_SHELL01 :: computeBmatrixAt: calling of this function is not allowed"); }
    virtual void computeNmatrixAt(GaussPoint *, FloatMatrix &)
    { _error("TR_SHELL01 :: computeNmatrixAt: calling of this function is not allowed"); }

    /// @todo In time delete
protected:
    virtual void computeGaussPoints()
    { _error("TR_SHELL01 :: computeGaussPoints: calling of this function is not allowed"); }
    virtual void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
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
