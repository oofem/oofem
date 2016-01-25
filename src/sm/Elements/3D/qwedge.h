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

#ifndef qwedge_h
#define qwedge_h

#include "Elements/structural3delement.h"
#include "ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_QWedge_Name "qwedge"

namespace oofem {
class FEI3dWedgeQuad;

/**
 * This class implements an Quadratic 3d  15 - node structural finite element.
 * Each node has 3 degrees of freedom.
 *
 * @author Mikael Öhman (among others)
 */
class QWedge : public Structural3DElement, public SPRNodalRecoveryModelInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface

{
protected:
    static FEI3dWedgeQuad interpolation;

public:
    QWedge(int, Domain *);
    virtual ~QWedge() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 9; }

    virtual Interface *giveInterface(InterfaceType);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QWedge_Name; }
    virtual const char *giveClassName() const { return "QWedge"; }

};
} // end namespace oofem
#endif
