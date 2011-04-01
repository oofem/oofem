/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.h,v 1.3 2003/04/23 14:22:15 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#ifndef tr21_2d_supg_h
#define tr21_2d_supg_h


#include "supgelement2.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "fei2dtrquad.h"
#include "fei2dtrlin.h"

#include "primaryfield.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "leplic.h"
#include "levelsetpcs.h"

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Class representing 2d triangular element  with quadratic velocity
 * and linear pressure aproximations for solving incompressible fluid problems
 * with SUPG solver
 *
 */

class TR21_2D_SUPG : public SUPGElement2, public LevelSetPCSElementInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI2dTrQuad velocityInterpolation;
    static FEI2dTrLin pressureInterpolation;
    IntArray pressureDofManArray;




public:
    // constructor
    TR21_2D_SUPG(int, Domain *);
    ~TR21_2D_SUPG();                        // destructor

    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models. Child classes should overload this function only
     * if they can be used together with nonlocal materil (where nonlocal averaging over
     * surronding volume is used).
     * @returns nonzero if successful; zero otherwise
     */
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */


    // definition
    const char *giveClassName() const { return "TR21_2D_SUPG"; }
    classType                giveClassID() const { return SUPGElementClass; }
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }

    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;
    virtual void           giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int            computeNumberOfDofs(EquationID ut);
    IRResultType           initializeFrom(InputRecord *ir);
    virtual void          updateYourself(TimeStep *tStep);
    /// used to check consistency and initialize some element geometry data (area,b,c)
    virtual int           checkConsistency();

    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType   saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered
     */
    contextIOResultType   restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);




    /**
     * @name The element interface required by LevelSetPCSElementInterface
     */

    /** Evaluetes F in level set equation of the form
     *  fi_t+F(grad(fi), x)*norm(grad(fi)) = 0
     *  where for interface position driven by flow with speed u:
     *  F=dotProduct(u,grad(fi))/norm(grad(fi))
     */
    virtual double LS_PCS_computeF(LevelSetPCS *, TimeStep *);

    /** Returns gradient of shape functions (assumed constatnt <- linear approx)
     */
    virtual void LS_PCS_computedN(FloatMatrix &answer);
    /// Returns receiver's volume
    virtual double LS_PCS_computeVolume();

    void LS_PCS_computeVolume(double &answer,  const FloatArray **coordinates);

    /** Evaluetes S in level set equation of the form
     *  fi_t = S(fi)*(1-norm(grad(fi))) = 0
     *  where for interface position driven by flow with speed u:
     *  S=fi/sqrt(fi^2+eps^2)
     */
    virtual double LS_PCS_computeS(LevelSetPCS *, TimeStep *);

    /**
     * Returns VOF fractions for each material on element
     * according to nodal values of level set function (passed as parameter)
     */
    virtual void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi);


    //@}
    /**
     * @name The element interface required by ZZNodalRecoveryModel
     */
    //@{
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * @param type determines the type of internal variable to be recovered
     * @return size of DofManger record required to hold recovered values
     */
    int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    /**
     * Returns the corresponding element to interface
     */
    Element *ZZNodalRecoveryMI_giveElement() { return this; }
    /**
     * Evaluates N matrix (interpolation estimated stress matrix).
     */
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint,
                                                             InternalStateType type);
    //@}


    /**
     * @name The element interface required by NodalAveragingRecoveryModel
     */
    //@{
    /**
     * Computes the element value in given node.
     * @param answer contains the result
     * @param node element node number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                    InternalStateType type, TimeStep *tStep);
    /**
     * Computes the element value in given side.
     * @param answer contains the result
     * @param node element side number
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                   InternalStateType type, TimeStep *tStep);
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * @param type determines the type of internal variable to be recovered
     * @return size of DofManger record required to hold recovered values
     */
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    //@}







    /*another helping functions for computing VOFFractions*/

    void computeIntersection(int iedge, FloatArray &intcoords, FloatArray &fi);

    void computeMiddlePointOnParabolicArc(FloatArray &answer, int iedge, FloatArray borderpoints);

    void computeCenterOf(FloatArray &C, FloatArray c, int dim);

    void computeQuadraticRoots(FloatArray Coeff, double &r1, double &r2);

    void computeCoordsOfEdge(FloatArray &answer, int iedge);

    void computeQuadraticFunct(FloatArray &answer, int iedge);

    void computeQuadraticFunct(FloatArray &answer, FloatArray line);

    /**
     * Returns the integration point corresponding value in REDUCED form.
     * @param answer contain corresponding ip value, zero sized if not available.
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    //
    // Graphics output
    //
    //void          drawYourself (oofegGraphicContext&);
    virtual void  drawRawGeometry(oofegGraphicContext &);
    virtual void  drawScalar(oofegGraphicContext &context);
    //virtual void  drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

    /** Prints output of receiver to stream, for given time step */
    virtual void   printOutputAt(FILE *, TimeStep *);
    double computeCriticalTimeStep(TimeStep *tStep);

    // three terms for computing their norms due to computing t_supg
    virtual void computeAdvectionTerm(FloatMatrix &answer, TimeStep *atTime);

    virtual void computeAdvectionDeltaTerm(FloatMatrix &answer, TimeStep *atTime);

    virtual void computeMassDeltaTerm(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeLSICTerm(FloatMatrix &answer, TimeStep *atTime);



    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);


protected:

    virtual void giveLocalVelocityDofMap (IntArray &map);
    virtual void giveLocalPressureDofMap (IntArray &map);

    void                  computeGaussPoints();
    virtual void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp);
    virtual void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual int  giveNumberOfSpatialDimensions();
    double  computeVolumeAround(GaussPoint *aGaussPoint);
    virtual void initGeometry();



    virtual void updateStabilizationCoeffs(TimeStep *);

    virtual int giveTermIntergationRuleIndex(CharType termType);
};
} // end namespace oofem
#endif // tr21_2d_supg_h
