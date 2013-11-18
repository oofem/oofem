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

//   ***************************************************************************************
//   *** 2D LINEAR TRIANGULAR ELEMENT FOR FLUID DYNAMIC PROBLEMS SOLVED WITH PFEM METHOD ***
//   ***************************************************************************************

#ifndef tr1_2d_pfem_h
#define tr1_2d_pfem_h


#include "pfemelement2d.h"
#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"




#include "fei2dtrlin.h"

///@name Input fields for TR1PFEM
//@{
#define _IFT_TR1_2D_PFEM_Name "tr1pfem"
//@}
namespace oofem {

class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This class is the implementation of triangular PFEM element with linear (and equal order) interpolation of velocity and pressure fields.   
 */
 class TR1_2D_PFEM : public PFEMElement2d
{
 protected:
  //double a[3];
  // ?????????????????????
  double b [ 3 ];
  double c [ 3 ];
  double area;


  
  static FEI2dTrLin velocityInterpolation;
  static FEI2dTrLin pressureInterpolation;

  static IntArray ordering;
  static IntArray edge_ordering [ 3 ];
public:
    /// Constructor - not used
    TR1_2D_PFEM(int, Domain *);
	/// Constructor - used
    TR1_2D_PFEM(int, Domain *, int, int, int, int, int);
	// Destructor
    ~TR1_2D_PFEM();

    virtual void          computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *);
    virtual void          computeDiagonalMassMtrx(FloatArray &answer, TimeStep *);
    virtual void          computeDiagonalMassMtrx(FloatMatrix &answer, TimeStep *);
    
    // NOT IN USE
    virtual void computeInvertedStabilizationDiagonalMassMtrx(FloatMatrix &answer, TimeStep *atTime);

    double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &inDirection);
    void computeStabilizationParameters (FloatArray &answer, GaussPoint *gp, TimeStep *atTime);

    /// calculates critical time step
    virtual double        computeCriticalTimeStep(TimeStep *tStep);

    /**
     * Computes the global coordinates from given element's local coordinates.
     * Required by nonlocal material models. Child classes should overload this function only
     * if they can be used together with nonlocal materil (where nonlocal averaging over
     * surronding volume is used).
     * @returns nonzero if successful; zero otherwise
     */
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    /**
     * Computes the element local coordinates from given global coordinates.
     * @returns nonzero if successful (if point is inside element); zero otherwise
     */
    virtual bool computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);

    // definition
    const char *giveClassName() const { return "PFEMElement"; }
    virtual const char *giveInputRecordName() const { return _IFT_TR1_2D_PFEM_Name; }

    classType                giveClassID() const { return PFEMElementClass; }
    Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;
	
    virtual void           giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int            computeNumberOfDofs();
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

    double computeVolumeAround(GaussPoint *aGaussPoint);

    /** Interface requesting service */
    Interface *giveInterface(InterfaceType);
    
    virtual Element *giveElement() { return this; }

    //@}

    /**
     * Returns the integration point corresponding value in REDUCED form.
     * @param answer contain corresponding ip value, zero sized if not available.
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    //virtual int giveIPValueSize(InternalStateType type, GaussPoint *);
    //virtual InternalStateValueType giveIPValueType(InternalStateType type);
    //virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);


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

    virtual FEInterpolation * giveVelocityInterpolation(){return & velocityInterpolation;}
    virtual FEInterpolation * givePressureInterpolation(){return & pressureInterpolation;}


 protected:
       
    void computeGaussPoints();
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *);
      
    virtual void computeForceVector(FloatArray &answer, TimeStep *atTime);//F

    // NOT IN USE
    virtual void computePFEMSubstitutionMatrix(FloatMatrix &answer, TimeStep *atTime );//S


	/// Calculates the body load vector
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *atTime);
	/// Calculates the boundary condition sub-vector on an edge
    virtual void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *atTime);

};

} // end namespace oofem
#endif // tr1_2d_pfem_h







