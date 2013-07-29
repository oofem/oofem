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

#ifndef qbrick1_ht_h
#define qbrick1_ht_h

#include "transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "eleminterpmapperinterface.h"
#include "fei3dhexaquad.h"

#define _IFT_QBrick1_ht_Name "qbrick1ht"
#define _IFT_QBrick1_hmt_Name "qbrick1hmt"

namespace oofem {
/**
 * Brick (3d) elements with quadratic approximation for heat and mass transfer. Each node has 1 (heat) or 2 (heat+moisture) degrees of freedom.
 * @author Vit Smilauer
 */
class QBrick1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface
{
protected:
    static FEI3dHexaQuad interpolation;
    int numberOfGaussPoints;

public:
    QBrick1_ht(int n, Domain *d);
    virtual ~ QBrick1_ht();

    virtual double computeVolumeAround(GaussPoint *gp);
    virtual FEInterpolation *giveInterpolation() const { return & interpolation; }
    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QBrick1_ht_Name; }
    virtual const char *giveClassName() const { return "QBrick1_ht"; }
    virtual classType giveClassID() const { return QBrick1_htClass; }

    virtual int computeNumberOfDofs(EquationID ut) { return ( emode == HeatTransferEM ) ? 20 : 40; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _3dHeat; }

    virtual Interface *giveInterface(InterfaceType t);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();
    
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

protected:
    virtual void computeGaussPoints();
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual IntegrationRule *GetSurfaceIntegrationRule(int approxOrder);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iEdge);
};

class QBrick1_hmt : public QBrick1_ht
{
public:
    QBrick1_hmt(int n, Domain *d);

    virtual MaterialMode giveMaterialMode() { return _3dHeMo; }
    virtual const char *giveInputRecordName() const { return _IFT_QBrick1_hmt_Name; }
    virtual const char *giveClassName() const { return "QBrick1_hmt"; }
    virtual classType giveClassID() const { return QBrick1_hmtClass; }
};

} // end namespace oofem
#endif // qbrick1_ht_h
