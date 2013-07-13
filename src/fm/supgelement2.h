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

#ifndef supgelement2_h
#define supgelement2_h

#include "supgelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"

#include "primaryfield.h"
#include "fluiddynamicmaterial.h"

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This abstract class represent a general base element class for
 * fluid dynamic problems. It is derived from SUPGElement base,
 * but it provides general algorithms how to obtain characteristic
 * vectors & matrices by means of numerical integration.
 * It is therefore suitable for high-order elements where exact
 * integration would be difficult.
 */
class SUPGElement2 : public SUPGElement
{
public:
    SUPGElement2(int n, Domain *aDomain);
    virtual ~SUPGElement2();

    virtual IRResultType initializeFrom(InputRecord *ir);

    // characteristic  matrix
    virtual void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *tStep);
    virtual double giveCharacteristicValue(CharType, TimeStep *tStep);
    virtual void updateElementForNewInterfacePosition(TimeStep *tStep) { }

    virtual void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime);
    virtual void computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime);
    virtual void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime);
    virtual void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime);
    virtual void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime);

    virtual double computeCriticalTimeStep(TimeStep *tStep) = 0;

    // time step termination
    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual int checkConsistency();

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    // Graphics output
    //void drawYourself(oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext&) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    /**
     * Returns index of integration rule used to evaluate given characteristic term.
     */
    virtual int giveTermIntergationRuleIndex(CharType termType) = 0;

    virtual void computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime) = 0;
    virtual void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp) = 0;
    virtual void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime) = 0;
    virtual void computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime) = 0;
    virtual int giveNumberOfSpatialDimensions() = 0;

    //virtual void computeEdgeNuMatrix (FloatMatrix& answer, GaussPoint* gp) = 0;
    //virtual double computeEdgeVolumeAround (GaussPoint*, int iedge) = 0;
    //virtual void computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iedge) = 0;
    //virtual int computeLoadGToLRotationMtrx (FloatMatrix& answer) = 0;
    //virtual int computeLoadLEToLRotationMatrix (FloatMatrix& answer, int iedge, GaussPoint* gp) = 0;
    //virtual void giveEdgeUDofMapping (IntArray& answer, int iedge) = 0;
    //virtual int giveUApproxOrder () = 0;
    virtual void computeEdgeLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *atTime);
    virtual void computeSurfaceLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *atTime);
    virtual void computeEdgeLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *atTime);
    virtual void computeSurfaceLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *atTime);
};
} // end namespace oofem
#endif // supgelement2_h
