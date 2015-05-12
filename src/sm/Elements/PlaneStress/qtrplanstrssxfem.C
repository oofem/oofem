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

#include "qtrplanstrssxfem.h"

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "vtkxmlexportmodule.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentitem.h"
#include "xfem/delaunay.h"
#include "xfem/XFEMDebugTools.h"
#include "feinterpol.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "domain.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "xfem/patchintegrationrule.h"
#endif



#include <string>
#include <sstream>

namespace oofem {
REGISTER_Element(QTrPlaneStress2dXFEM);

QTrPlaneStress2dXFEM::~QTrPlaneStress2dXFEM() {
	// TODO Auto-generated destructor stub
}

void QTrPlaneStress2dXFEM :: updateYourself(TimeStep *tStep)
{
    QTrPlaneStress2d :: updateYourself(tStep);
    XfemStructuralElementInterface :: updateYourselfCZ(tStep);
}

void QTrPlaneStress2dXFEM :: postInitialize()
{
    QTrPlaneStress2d :: postInitialize();
    XfemStructuralElementInterface :: initializeCZMaterial();
}

int QTrPlaneStress2dXFEM :: computeNumberOfDofs()
{
    int ndofs = 0;

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        ndofs += this->giveDofManager(i)->giveNumberOfDofs();
    }

    return ndofs;
}

void QTrPlaneStress2dXFEM :: computeGaussPoints()
{
    if ( giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = giveDomain()->giveXfemManager();

        if ( xMan->isElementEnriched(this) ) {
            if ( !this->XfemElementInterface_updateIntegrationRule() ) {
                QTrPlaneStress2d :: computeGaussPoints();
            }
        } else {
            QTrPlaneStress2d :: computeGaussPoints();
        }
    } else   {
        QTrPlaneStress2d :: computeGaussPoints();
    }
}

void QTrPlaneStress2dXFEM :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    XfemElementInterface_createEnrNmatrixAt(answer, iLocCoord, * this, false);
}

void QTrPlaneStress2dXFEM :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    XfemElementInterface_createEnrBmatrixAt(answer, * gp, * this);
}

void QTrPlaneStress2dXFEM :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    XfemElementInterface_createEnrBHmatrixAt(answer, * gp, * this);
}

void
QTrPlaneStress2dXFEM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // Continuous part
    QTrPlaneStress2d :: giveDofManDofIDMask(inode, answer);

    // Discontinuous part
    if( this->giveDomain()->hasXfemManager() ) {
        DofManager *dMan = giveDofManager(inode);
        XfemManager *xMan = giveDomain()->giveXfemManager();

        int placeInArray = domain->giveDofManPlaceInArray(dMan->giveGlobalNumber());
        const std::vector<int> &nodeEiIndices = xMan->giveNodeEnrichmentItemIndices( placeInArray );
        for ( size_t i = 0; i < nodeEiIndices.size(); i++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(nodeEiIndices[i]);
            if ( ei->isDofManEnriched(* dMan) ) {
                IntArray eiDofIdArray;
                ei->computeEnrichedDofManDofIdArray(eiDofIdArray, *dMan);
                answer.followedBy(eiDofIdArray);
            }
        }
    }
}

void
QTrPlaneStress2dXFEM :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface :: XfemElementInterface_computeConstitutiveMatrixAt(answer, rMode, gp, tStep);
}

void
QTrPlaneStress2dXFEM :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface :: XfemElementInterface_computeStressVector(answer, strain, gp, tStep);
}

void QTrPlaneStress2dXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    QTrPlaneStress2d :: computeStiffnessMatrix(answer, rMode, tStep);
    XfemStructuralElementInterface :: computeCohesiveTangent(answer, tStep);

    const double tol = 1.0e-6;
    const double regularizationCoeff = 1.0e-6;
    int numRows = answer.giveNumberOfRows();
    for(int i = 0; i < numRows; i++) {
        if( fabs(answer(i,i)) < tol ) {
            answer(i,i) += regularizationCoeff;
//          printf("Found zero on diagonal.\n");
        }
    }
}

void QTrPlaneStress2dXFEM :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    XfemStructuralElementInterface::XfemElementInterface_computeDeformationGradientVector(answer, gp, tStep);
}

void
QTrPlaneStress2dXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    QTrPlaneStress2d :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    XfemStructuralElementInterface :: computeCohesiveForces(answer, tStep);
}

Element_Geometry_Type
QTrPlaneStress2dXFEM :: giveGeometryType() const
{
    if ( this->giveDomain()->hasXfemManager() ) {
        XfemManager *xMan = this->giveDomain()->giveXfemManager();
        if ( xMan->isElementEnriched(this) ) {
            return EGT_Composite;
        } else {
            return EGT_Composite;
        }
    } else   {
        return EGT_triangle_2;
    }
}

IRResultType
QTrPlaneStress2dXFEM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    result = QTrPlaneStress2d :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    result = XfemStructuralElementInterface :: initializeCZFrom(ir);
    return result;
}

MaterialMode QTrPlaneStress2dXFEM :: giveMaterialMode()
{
    return XfemStructuralElementInterface :: giveMaterialMode();
}

void QTrPlaneStress2dXFEM :: giveInputRecord(DynamicInputRecord &input)
{
    QTrPlaneStress2d :: giveInputRecord(input);
    XfemStructuralElementInterface :: giveCZInputRecord(input);
}

void
QTrPlaneStress2dXFEM :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    // TODO: Validate implementation.

    FloatArray u;
    FloatMatrix n;

    XfemElementInterface_createEnrNmatrixAt(n, lcoords, * this, false);

    this->computeVectorOf(mode, tStep, u);
    answer.beProductOf(n, u);
}


} /* namespace OOFEM */
