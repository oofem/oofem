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

#ifndef tetrah1_ht_h
#define tetrah1_ht_h

#include "transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"

#define _IFT_Tetrah1_ht_Name "tetrah1ht"
#define _IFT_Tetrah1_hmt_Name "tetrah1hmt"

namespace oofem {
class FEI3dTetLin;

/**
 * Tetrahedral (3d) element with linear approximation for heat and mass transfer.
 */
class Tetrah1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI3dTetLin interpolation;

public:
    Tetrah1_ht(int n, Domain * d);
    virtual ~Tetrah1_ht();

    virtual FEInterpolation *giveInterpolation() const;

    virtual double computeVolumeAround(GaussPoint *gp);

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_Tetrah1_ht_Name; }
    virtual const char *giveClassName() const { return "Tetrah1_ht"; }

    virtual int computeNumberOfDofs() { return ( emode == HeatTransferEM ) ? 4 : 8; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _3dHeat; }

    virtual Interface *giveInterface(InterfaceType t);
    virtual int testElementExtension(ElementExtension ext)
    { return ( ( ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ) ) ? 1 : 0 ); }

#ifdef __OOFEG
    // Graphics output
    virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep);
    virtual void drawScalar(oofegGraphicContext &gc, TimeStep *tStep);
    //virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) {}
#endif

protected:
    virtual void computeGaussPoints();
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual IntegrationRule *GetSurfaceIntegrationRule(int approxOrder);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iEdge);
};

class Tetrah1_hmt : public Tetrah1_ht
{
public:
    Tetrah1_hmt(int n, Domain * d);

    virtual const char *giveInputRecordName() const { return _IFT_Tetrah1_hmt_Name; }
    virtual const char *giveClassName() const { return "Tetrah1_hmt"; }
    virtual MaterialMode giveMaterialMode() { return _3dHeMo; }
};
} // end namespace oofem
#endif // tetrah1_ht_h
