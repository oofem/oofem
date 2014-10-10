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

#ifndef quad1_ht_h
#define quad1_ht_h

#include "transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"

#define _IFT_Quad1_ht_Name "quad1ht"
#define _IFT_Quad1_hmt_Name "quad1hmt"
#define _IFT_Quad1_mt_Name "quad1mt"

namespace oofem {
class FEI2dQuadLin;

/**
 * Quadratic (2d) element with linear approximation for heat transfer.
 */
class Quad1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dQuadLin interpolation;

public:
    Quad1_ht(int n, Domain * d);
    virtual ~Quad1_ht();

    virtual FEInterpolation *giveInterpolation() const;
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual const char *giveInputRecordName() const { return _IFT_Quad1_ht_Name; }
    virtual const char *giveClassName() const { return "Quad1_ht"; }

    //    virtual int computeNumberOfDofs() { return ( emode == HeatTransferEM ) ? 4 : 8; }
    virtual int computeNumberOfDofs() { return 4; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _2dHeat; }
    virtual double giveThicknessAt(const FloatArray &gcoords);

    virtual Interface *giveInterface(InterfaceType t);

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
};

/**
 * Class for heat and mass transfer.
 */
class Quad1_hmt : public Quad1_ht
{
public:
    Quad1_hmt(int n, Domain * d);

    virtual const char *giveInputRecordName() const { return _IFT_Quad1_hmt_Name; }
    virtual const char *giveClassName() const { return "Quad1_hmt"; }
    virtual int computeNumberOfDofs() { return 8; }
    virtual MaterialMode giveMaterialMode() { return _2dHeMo; }
};

/**
 * Class for mass transfer.
 */
class Quad1_mt : public Quad1_ht
{
public:
    Quad1_mt(int n, Domain * d);

    virtual const char *giveInputRecordName() const { return _IFT_Quad1_mt_Name; }
    virtual const char *giveClassName() const { return "Quad1_mt"; }
    virtual int computeNumberOfDofs() { return 4; }
    virtual MaterialMode giveMaterialMode() { return _2dHeat; }
};
} // end namespace oofem
#endif // quad1_ht_h
