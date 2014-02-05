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

#ifndef brouzoulisshell_h
#define brouzoulisshell_h

#include "nlstructuralelement.h"
#include "fei2dquadlin.h"

#include "floatmatrix.h"
#include "floatarray.h"


#define _IFT_BrouzoulisShell_Name "BrouzoulisShell"

namespace oofem {



/**
 * This class implements an bi-linear quad shell element - (4 nodes) Each node has 5 or 6? degrees of freedom.
 * Based on th paper: Advances in one-point quadrature shell elements
 * By: Ted Belytschko, Bak Leong Wong and Huai-Yang Chiang
 * Computer Methods in Applied Mechanics and Engineering 96 (1992) 93-107
 *
 * @author Jim Brouzoulis
 */
class BrouzoulisShell : public NLStructuralElement
{

protected:
    static FEI2dQuadLin interpolation;

    bool matRotation;

    /// Last equilibrium triad at the centre (xi1 = xi2 = 0).
    FloatMatrix corotTriad;
    /// Temporary triad at the centre.
    FloatMatrix tempCorotTriad;
    
public:
    BrouzoulisShell(int n, Domain *d);
    virtual ~BrouzoulisShell() {}

    virtual FEInterpolation *giveInterpolation() const { return & interpolation; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);


    virtual Interface *giveInterface(InterfaceType);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_BrouzoulisShell_Name; }
    virtual const char *giveClassName() const { return "BrouzoulisShell"; }
    virtual int computeNumberOfDofs() { return 20; } // 4*5 dofs
    virtual int computeNumberOfGlobalDofs() { return 24; } // 4*6 dofs
    virtual MaterialMode giveMaterialMode() { return _2dPlate; };
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_quad_1; };


    // Shell specific
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void computeSectionalForcesAt(FloatArray &sectionalForces, FloatArray &stress, GaussPoint *gp, TimeStep *tStep);
    void computedNdX(FloatMatrix &answer);

    void setupInitialNodeDirectors();
    void giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    FloatArray &giveInitialNodeDirector(int i) {
        return this->initialNodeDirectors [ i - 1 ];
    }

    FloatArray &giveInitialNodeMidSurface(int i) {
        return this->initialNodeMidSurface [ i - 1 ];
    }


    void computeCoRotatedNodeCoords(FloatMatrix &answer);

protected:
    std :: vector< FloatArray >initialNodeDirectors;
    std :: vector< FloatArray >initialNodeMidSurface;
    virtual void computeGaussPoints();

    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual int giveApproxOrder() { return 2; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 35; }

    void createInternalNodes();
    virtual void postInitialize(){ 
        NLStructuralElement :: postInitialize();
        //createInternalNodes(); 
        //this->setupInitialNodeDirectors();
    };
    /**
     * @name Surface load support
     */
    //@{
    virtual IntegrationRule *GetSurfaceIntegrationRule(int);
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void giveSurfaceDofMapping(IntArray &answer, int) const;
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int);
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    //@}
};
} // end namespace oofem

#if 0
class BrouzoulisShell : public NLStructuralElement, NodalAveragingRecoveryModelInterface
{

protected:
    static FEI3dWedgeLin interpolationUV;
    static FEI3dWedgeQuad interpolationW;

    bool matRotation;

public:
    BrouzoulisShell(int n, Domain *d);
    virtual ~BrouzoulisShell() {}

    virtual FEInterpolation *giveInterpolation() const { return & interpolationUV; } ///@todo fix

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);


    virtual Interface *giveInterface(InterfaceType);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_BrouzoulisShell_Name; }
    virtual const char *giveClassName() const { return "BrouzoulisShell"; }
    virtual int computeNumberOfDofs() { return 27; } // 6*3 dofs + 9*1 dofs
    virtual MaterialMode giveMaterialMode() { return _3dMat; };
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_wedge_2; }

    // Nodal averaging interface:
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);


    // Shell specific
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void computeSectionalForcesAt(FloatArray &sectionalForces, FloatArray &stress, GaussPoint *gp, TimeStep *tStep);
    void computeLambdaMatrices(FloatMatrix lambda [ 3 ], double zeta);
    void evalInitialCovarBaseVectorsAt(FloatArray &lcoords, FloatMatrix &Gcov);
    void giveBmat(FloatArray &B, GaussPoint *gp);
    void setupInitialNodeDirectors();
    void giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords);

    FloatArray &giveInitialNodeDirector(int i) {
        return this->initialNodeDirectors [ i - 1 ];
    }

    FloatArray &giveInitialNodeMidSurface(int i) {
        return this->initialNodeMidSurface [ i - 1 ];
    }
protected:
    std :: vector< FloatArray >initialNodeDirectors;
    std :: vector< FloatArray >initialNodeMidSurface;
    virtual void computeGaussPoints();

    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    virtual void computeBHmatrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual int giveApproxOrder() { return 2; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 35; }

    void createInternalNodes();
    virtual void postInitialize(){ 
        NLStructuralElement :: postInitialize();
        //createInternalNodes(); 
        this->setupInitialNodeDirectors();
    };
    /**
     * @name Surface load support
     */
    //@{
    virtual IntegrationRule *GetSurfaceIntegrationRule(int);
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void giveSurfaceDofMapping(IntArray &answer, int) const;
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int);
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    //@}
};
#endif

#endif
