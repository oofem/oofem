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

#ifndef QTRPLANSTRSSXFEM_H_
#define QTRPLANSTRSSXFEM_H_

#include "../sm/Elements/PlaneStress/qtrplstr.h"
#include "../sm/xfem/xfemstructuralelementinterface.h"
#include "vtkxmlexportmodule.h"

#include "domain.h"

#define _IFT_QTrPlaneStress2dXFEM_Name "qtrplanestress2dxfem"

namespace oofem {


/**
 * 6-node triangle with XFEM kinematics.
 * Only corner nodes are enriched, i.e. we use a quadratic approximation
 * for the continuous field and a linear approximation for enrichment fields.
 * @author Erik Svenning
 * @date Dec 19, 2014
 */
class QTrPlaneStress2dXFEM : public QTrPlaneStress2d, public XfemStructuralElementInterface, public VTKXMLExportModuleElementInterface
{
protected:
    virtual void updateYourself(TimeStep *tStep);
    virtual void postInitialize();

public:
    QTrPlaneStress2dXFEM(int n, Domain * d) : QTrPlaneStress2d(n, d), XfemStructuralElementInterface(this), VTKXMLExportModuleElementInterface() { numberOfDofMans = 6; }
    virtual ~QTrPlaneStress2dXFEM();

    virtual const char *giveInputRecordName() const { return _IFT_QTrPlaneStress2dXFEM_Name; }
    virtual const char *giveClassName() const { return "QTrPlaneStress2dXFEM"; }


    virtual int computeNumberOfDofs();
    virtual void computeGaussPoints();
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                  int lowerIndx = 1, int upperIndx = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL) { XfemStructuralElementInterface :: XfemElementInterface_computeConsistentMassMatrix(answer, tStep, mass, ipDensity); }

    virtual Element_Geometry_Type giveGeometryType() const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode();
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer);

protected:
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 6; } // TODO: Check

};

} /* namespace OOFEM */
#endif /* QTRPLANSTRSSXFEM_H_ */
