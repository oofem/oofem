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

#ifndef tr1_ht_h
#define tr1_ht_h

#include "transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"

#define _IFT_Tr1_hmt_Name "tr1hmt"
#define _IFT_Tr1_ht_Name "tr1ht"


namespace oofem {
class FEI2dTrLin;

/**
 * Triangle (2d) element with linear approximation for heat transfer.
 * @todo Use the interpolation classes.
 */
class Tr1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dTrLin interp;

public:
    Tr1_ht(int n, Domain * d);
    virtual ~Tr1_ht();

    virtual double computeVolumeAround(GaussPoint *gp);

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_Tr1_ht_Name; }
    virtual const char *giveClassName() const { return "Tr1_htElement"; }

    virtual int computeNumberOfDofs() { return ( emode == HeatTransferEM ) ? 3 : 6; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _2dHeat; }
    virtual double giveThicknessAt(const FloatArray &gcoords);

    virtual Interface *giveInterface(InterfaceType t);

    virtual FEInterpolation *giveInterpolation() const;

#ifdef __OOFEG
    // Graphics output
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
class Tr1_hmt : public Tr1_ht
{
public:
    Tr1_hmt(int n, Domain * d);

    virtual const char *giveInputRecordName() const { return _IFT_Tr1_hmt_Name; }
    virtual const char *giveClassName() const { return "Tr1_hmt"; }
    virtual MaterialMode giveMaterialMode() { return _2dHeMo; }
};
} // end namespace oofem
#endif // tr1_ht_h
