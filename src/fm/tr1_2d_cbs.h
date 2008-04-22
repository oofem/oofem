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

//   *****************************************************************************************
//   *** 2D LINEAR TRIANGULAR ELEMENT FOR FLUID DYNAMIC PROBLEMS SOLVED WITH CBS ALGORITHM ***
//   *****************************************************************************************

#ifndef tr1_2d_cbs_h
#define tr1_2d_cbs_h


#include "cbselement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "primaryfield.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
//<RESTRICTED_SECTION>
#include "leplic.h"
//</RESTRICTED_SECTION>

class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This class is the implementation of triangular CFD element with linear (and equal order) interpolation of velocity and pressure fields.
 * Should be used with CBS solution algorithm.
 */
class TR1_2D_CBS : public CBSElement, public SpatialLocalizerInterface, public EIPrimaryFieldInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface, public SPRNodalRecoveryModelInterface
    //<RESTRICTED_SECTION>
    , public LEPlicElementInterface
    //</RESTRICTED_SECTION>
{
protected:
    //double a[3];
    double b [ 3 ];
    double c [ 3 ];
    double area;

public:
    // constructor
    TR1_2D_CBS(int, Domain *);
    ~TR1_2D_CBS();                        // destructor

    /** Calculates consistent mass matrix */
    virtual void          computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *);
    /** Calculates diagonal mass matrix */
    virtual void          computeDiagonalMassMtrx(FloatArray &answer, TimeStep *);
    /** Calculates convection component for (*) velocities */
    virtual void          computeConvectionTermsI(FloatArray &answer, TimeStep *);
    /** Calculates contribution from diffusion terms for (*) velocities */
    virtual void          computeDiffusionTermsI(FloatArray &answer, TimeStep *);
    /// computes velocity terms on RHS for density equation
    virtual void          computeDensityRhsVelocityTerms(FloatArray &answer, TimeStep *tStep);
    /// computes pressure terms on RHS for density equation
    virtual void          computeDensityRhsPressureTerms(FloatArray &answer, TimeStep *tStep);
    /// computes prescribed pressure due to applied tractions
    virtual void          computePrescribedTractionPressure(FloatArray &answer, TimeStep *tStep);
    /// computes number of edges/sides with prescribed traction contributing to node with prescribed pressure
    virtual void          computeNumberOfNodalPrescribedTractionPressureContributions(FloatArray &answer, TimeStep *tStep);
    /// calculates the pressure LHS
    virtual void          computePressureLhs(FloatMatrix &answer, TimeStep *tStep);
    /// calculates the RHS of velocity correction step
    virtual void          computeCorrectionRhs(FloatArray &answer, TimeStep *tStep);
    /// calculates critical time step
    virtual double        computeCriticalTimeStep(TimeStep *tStep);

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
    const char *giveClassName() const { return "CBSElement"; }
    classType                giveClassID() const { return CBSElementClass; }
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
    virtual void EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                       FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                       TimeStep *atTime);
    //@}

    //<RESTRICTED_SECTION>
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
    virtual void giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool upd);
    /// Assembles polygon representing receiver
    virtual void formMyVolumePoly(Polygon &myPoly, LEPlic *mat_interface, bool updFlag);
    virtual Element *giveElement() { return this; }
    virtual double computeMyVolume(LEPlic *matInterface, bool updFlag);
    /// Critical time step; limits already handled by computeCriticalTimeStep
    virtual double computeCriticalLEPlicTimeStep(TimeStep *tStep) { return 1.e6; }

    //@}
    //</RESTRICTED_SECTION>
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
    void                  computeGaussPoints();
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *);
};

#endif // tr1_2d_cbs_h







