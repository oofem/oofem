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


#ifndef lwedge_h
#define lwedge_h

#include "nlstructuralelement.h"
#include "fei3dwedgelin.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_LWedge_Name "lwedge"

namespace oofem {

/**
 * This class implements a Linear 3d  6 - node elasticity finite element.
 * 
 * Each node has 3 degrees of freedom.
 * One single additional attribute is needed for Gauss integration purpose :
 * 'jacobianMatrix'. This 3x3 matrix contains polynomials.
 */    
class LWedge : public NLStructuralElement, public SPRNodalRecoveryModelInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    int numberOfGaussPoints;
    static FEI3dWedgeLin interpolation;

public:
    LWedge(int, Domain *);
    virtual ~LWedge() {}

    virtual FEInterpolation *giveInterpolation() const { return &interpolation; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *);
   
    virtual Interface *giveInterface(InterfaceType);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }
    virtual int giveApproxOrder() { return 1; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 9; }


    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);

    //
    // definition & identification
    //
    virtual const char *giveInputRecordName() const { return _IFT_LWedge_Name; }
    virtual const char *giveClassName() const { return "LWedge"; }
    virtual classType giveClassID() const { return LWedgeClass; }
    virtual int computeNumberOfDofs(EquationID ut) { return 18; }
    virtual MaterialMode giveMaterialMode();

protected:
    virtual void computeGaussPoints();
    virtual void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *, FloatMatrix &);
};
} // end namespace oofem
#endif
