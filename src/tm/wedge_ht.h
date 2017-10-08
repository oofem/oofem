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
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of.
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#ifndef wedge_ht_h
#define wedge_ht_h


#include "transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"


#define _IFT_Wedge_ht_Name "wedgeht"
#define _IFT_Wedge_hmt_Name "wedgehmt"
#define _IFT_Wedge_mt_Name "wedgemt"

namespace oofem {
class FEI3dWedgeLin;

/**
 * This class implements a Linear 3d  6 - node thermal finite element.
 */
 class Wedge_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface, public SPRNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI3dWedgeLin interpolation;

public:
    Wedge_ht(int, Domain *);
    virtual ~Wedge_ht() { }


    virtual double computeVolumeAround(GaussPoint *gp);
    virtual FEInterpolation *giveInterpolation() const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Wedge_ht_Name; }
    virtual const char *giveClassName() const { return "Wedge_ht"; }

    virtual int computeNumberOfDofs() { return 6; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _3dHeat; }

    virtual Interface *giveInterface(InterfaceType t);
    virtual int testElementExtension(ElementExtension ext)
    { return ( ext == Element_EdgeLoadSupport );}

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);

 protected:
    virtual void computeGaussPoints();
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);

    
};



/**
 * Class for heat and mass transfer.
 */
class Wedge_hmt : public Wedge_ht
{
public:
    Wedge_hmt(int n, Domain * d);

    virtual const char *giveInputRecordName() const { return _IFT_Wedge_hmt_Name; }
    virtual const char *giveClassName() const { return "Wedge_hmt"; }
    virtual int computeNumberOfDofs() { return 12; }
    virtual MaterialMode giveMaterialMode() { return _3dHeMo; }
};

/**
 * Class for mass transfer.
 */
class Wedge_mt : public Wedge_ht
{
public:
    Wedge_mt(int n, Domain * d);

    virtual const char *giveInputRecordName() const { return _IFT_Wedge_mt_Name; }
    virtual const char *giveClassName() const { return "Wedge_mt"; }
    virtual int computeNumberOfDofs() { return 6; }
    virtual MaterialMode giveMaterialMode() { return _3dHeat; }
};

 
} // end namespace oofem
#endif
