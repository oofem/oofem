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

//   ************************************************************************************************
//   *** 2D LINEAR TRIANGULAR ELEMENT FOR FLUID DYNAMIC PROBLEMS SOLVED WITH SUPG/PSSPG ALGORITHM ***
//   ************************************************************************************************

#ifndef tr1_2d_supg_h
#define tr1_2d_supg_h


#include "supgelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"

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
 * Class representing 2d linear triangular element
 * for solving incompressible fluid with SUPG solver
 *
 * This class is similar to TR1_2D_SUPG2, but diference is in handling
 * multiple fluids. This class uses rule of mixture which interpolates the
 * properties using VOF value, requiring the use of twofluidmaterial class
 * as material model for this situation.
 */
class TR1_2D_SUPG : public SUPGElement, public SpatialLocalizerInterface, public EIPrimaryFieldInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface
    , public LEPlicElementInterface, public LevelSetPCSElementInterface
{
protected:
    //double a[3];
    double b [ 3 ];
    double c [ 3 ];
    double area;

    /// Cached transformation matrix of receiver
    FloatMatrix rotationMatrix;
    /// Flag indicating if tranformation matrix has been already computed
    int rotationMatrixDefined;

public:
    // constructor
    TR1_2D_SUPG(int, Domain *);
    ~TR1_2D_SUPG();                        // destructor

    /**
     * Computes acceleration terms (generalized mass matrix with stabilization terms ) for momentum balance equations(s)
     */
    void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes nonlinear advection terms for momentum balance equations(s)
     */
    void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime);
    /**
     * Computes the derivative of advection terms for momentum balance equations(s)
     * with respect to nodal velocities
     */
    void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /**
     *  Computes diffusion terms for momentum balance equations(s)
     */
    void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime);
    /** Computes the derivative of diffusion terms for momentum balance equations(s)
     *  with respect to nodal velocities
     */
    void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime);
    /** Computes pressure terms for momentum balance equations(s) */
    void computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /** Computes SLIC stabilization term for momentum balance equation(s) */
    void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /** Computes the linear advection term for mass conservation equation */
    void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes advection terms for mass conservation equation
     */
    void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime);
    /** Computes the derivative of advection terms for mass conservation equation
     *  with respect to nodal velocities
     */
    void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes diffusion terms for mass conservation equation
     */
    void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime) {
        answer.resize(3, 6);
        answer.zero();
    }
    void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime) {
        answer.resize(3);
        answer.zero();
    }
    /**
     * Computes acceleration terms for mass conservation equation
     */
    void  computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes pressure terms for mass conservation equation
     */
    void computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    // calculates critical time step
    // virtual double        computeCriticalTimeStep (TimeStep* tStep);
    /**
     * Computes Rhs terms due to boundary conditions
     */
    void  computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime);
    /**
     * Computes Rhs terms due to boundary conditions
     */
    void  computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime);


    void     updateStabilizationCoeffs(TimeStep *);
    /// calculates critical time step
    double        computeCriticalTimeStep(TimeStep *tStep);

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
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    // definition
    const char *giveClassName() const { return "SUPGElement"; }
    classType                giveClassID() const { return SUPGElementClass; }
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

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


    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);

    /**
     * @name The element interface required by SpatialLocalizerInterface
     */
    //@{
    /// Returns reference to corresponding element
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    /// Returns nonzero if given element contains given point
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    /// Returns distance of given point from element parametric center
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);
    //@}

    /**
     * @name The element interface required by SpatialLocalizerInterface
     */
    //@{
    /**
     * Evaluates the value of field at given point of interest (should be located inside receiver's volume) using
     * element interpolation.
     */
    virtual int  EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                       FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                       TimeStep *atTime);
    //@}

    /**
     * @name The element interface required by LEPlicElementInterface
     */
    //@{
    /// Computes corresponding volume fraction to given interface position
    virtual double computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag);
    /// Assembles the element material polygon
    virtual void formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                        const FloatArray &normal, const double p, bool updFlag);
    /// Assembles receiver material polygon based solely on given interface line
    virtual void formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                         const FloatArray &normal, const double p, bool updFlag);

    /// Truncates given material polygon to receiver
    virtual double truncateMatVolume(const Polygon &matvolpoly, double &volume);
    /// Computes the receiver center (in updated Lagrangian configuration)
    virtual void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool updFlag);
    /// Assembles polygon representing receiver
    virtual void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag);
    virtual Element *giveElement() { return this; }
    virtual double computeMyVolume(LEPlic *matInterface, bool updFlag);
    virtual double computeCriticalLEPlicTimeStep(TimeStep *tStep);

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


    /**
     * @name The element interface required by SPRNodalRecoveryModelInterface
     */
    //@{
    void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
    { return ZZNodalRecoveryMI_giveDofManRecordSize(type); }
    int SPRNodalRecoveryMI_giveNumberOfIP();
    //void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type);
    void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    SPRPatchType SPRNodalRecoveryMI_givePatchType();
    //@}

    /**
     * @name The element interface required by LevelSetPCSElementInterface
     */
    //@{
    /** Evaluetes F in level set equation of the form
     *  fi_t+F(grad(fi), x)*norm(grad(fi)) = 0
     *  where for interface position driven by flow with speed u:
     *  F=dotProduct(u,grad(fi))/norm(grad(fi))
     */
    double LS_PCS_computeF(LevelSetPCS *, TimeStep *);

    /** Returns gradient of shape functions (assumed constatnt <- linear approx)
     */
    void LS_PCS_computedN(FloatMatrix &answer);
    ///Returns receiver's volume
    double LS_PCS_computeVolume() { return area; }
    virtual double LS_PCS_computeS(LevelSetPCS *, TimeStep *);
    void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi);
    //@}

    /**
     * Transformation matrices
     * Updates rotation matrix r(l)=T r(g*) between  local and global coordinate system
     * taking into account also possible local - coordinate system in some elements
     * nodes.
     * Default implementation uses \ref computeGtoLRotationMatrix and
     * \ref computeGNDofRotationMatrix  services to compute result.
     * Default implembalentation uses cached rotation matrix in
     * rotationMatrix attribute, so rotation matrix is computed only once.
     * @return nonzero if transformation is necessary.
     */
    virtual int   updateRotationMatrix();      //
    /**
     * Returns  transformation matrix from global coord. system to local element
     * coordinate system ( i.e. r(l)=T r(g)). If no trasformation is necessary
     * then answer is empty mtrx and zero value is returned.
     * @return nonzero if transformation is necessary, zero otherwise.
     */
    virtual int  computeGtoLRotationMatrix(FloatMatrix &answer) {
        answer.beEmptyMtrx();
        return 0;
    }
    /**
     * Returns transformation matrix for DOFs from global coordinate system
     * to local coordinate system in nodes (i.e. r(n)=T r(g)) if mode == _toNodalCS.
     * If mode == _toGlobalCS, the transformation from local nodal cs to
     * global cs in node is returned. If no trasformation is
     * necessary sets answer to empty mtrx and returns zero value.
     * @return nonzero if transformation is necessary, zero otherwise.
     */
    virtual int  computeGNDofRotationMatrix(FloatMatrix &answer, DofManTransfType mode);
    virtual int computeGNLoadRotationMatrix(FloatMatrix &answer, DofManTransfType mode);

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


protected:
    virtual void giveLocalVelocityDofMap(IntArray &map);
    virtual void giveLocalPressureDofMap(IntArray &map);

    void                  computeNMtrx(FloatArray &answer, GaussPoint *gp);
    void                  computeGaussPoints();
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *);
    virtual void initGeometry();
};
} // end namespace oofem
#endif // tr1_2d_supg_h
