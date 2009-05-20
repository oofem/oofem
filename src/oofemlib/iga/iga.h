/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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

#ifndef iga_h
#define iga_h

/*
  oofem nodes - controll points (coordinates ) + dofs
  oofem elements - NURBS patches as well as integration elements
  
  
  NURBS PATCH:
  knot vector - store knot coordinates + multiplicity
  patch integration rule - keep list of elements
  
  integration element
  
  FEInterpolation:
  - need to be enriched, as one should (or can) pass knot span to evaluation routines 
  (in this way the patch by patch evaluation can be faster)
  
*/

#include "inputrecord.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "structuralcrosssection.h"
#include "mathfem.h"

class BSplineInterpolation : public FEInterpolation  
{
protected:
	/// number of spatial directions
	int nsd;        
	/// degree in each direction 
	int *degree;                               // eg. 2
	/// numberOfControllPoints[nsd]
	int *numberOfControllPoints; 
	/// knotVectors[nsd][sd_controll_points]
	double **knotVector;                       // eg. 0 0 0 1 2 3 4 4 5 5 5
	/// nonzero spans in each directions[nsd]
	int *numberOfKnotSpans;                    // eg. 5 (0-1,1-2,2-3,3-4,4-5)
	/// knot multiplicity
	int ** knotMultiplicity;                   // eg. 3 1 1 1 2 3
public:
  BSplineInterpolation (int nsd) : FEInterpolation (0) {this->nsd=nsd;}
  ~BSplineInterpolation();

  IRResultType initializeFrom(InputRecord *ir);                   
  virtual int giveNumberOfKnotSpans(int dim) {return numberOfKnotSpans[dim-1];}
  virtual double ** const  giveKnotVector() {return this->knotVector;}
  virtual int** const giveKnotMultiplicity() {return this->knotMultiplicity;}
  /**
   * Evaluates the array of interpolation functions (shape functions) at given point.
   * @param answer contains resulting array of evaluated interpolation functions
   * @param lcoords array containing (local) coordinates
   * @param span identifies local sub-patch id (span coordinates)
   * @param time time
   *
   * see also giveNonzeroBasisFunctMask method
   */
  virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const IntArray *span,  double time) ;
  virtual void evalN(FloatArray &answer, const FloatArray &lcoords, double time) ;
	/**
	 * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
	 * These derivatives are in global coordinate system (where the nodal coordinates are defined)
	 * @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
	 * @param coords coordinates of nodes defining the interpolation geometry
	 * @param lcoords array containing (local) coordinates
	 * @param time time
	 */
  virtual void evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, const IntArray* knotSpan, double time) ;
  virtual void evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, double time) ;
	virtual void evaldNdx(FloatMatrix &answer, Domain *d, IntArray &nodes, const FloatArray &lcoords, const IntArray* span, double time) ;
	virtual void evaldNdx(FloatMatrix &answer, Domain *d, IntArray &nodes, const FloatArray &lcoords, double time) ;
  /**
   * Evaluates global coordinates from given local ones
   * @param answer contains resulting global coordinates
   * @param nodes array of node numbers defining the interpolation geometry
   * @param lcoords array containing (local) coordinates
   * @param time time
   */
  virtual void local2global(FloatArray &answer, Domain *d, IntArray &nodes, const FloatArray &lcoords, double time) ;
  virtual void local2global(FloatArray &answer, const FloatArray**, const FloatArray &lcoords, double time) {
    OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
  }
  virtual int  global2local(FloatArray &answer, Domain *d, IntArray &nodes, const FloatArray &lcoords, double time) {
    OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
    return 0;
  }
  virtual int  global2local(FloatArray &answer, const FloatArray **coords, const FloatArray &lcoords, double time) {
    OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
    return 0;
  }
  /**
   * Evaluates the jacobian of transformation between local and global coordinates.
   */
  virtual double giveTransformationJacobian(const FloatArray **coords, const FloatArray &lcoords, double time) {
    OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
    return 0.0;
  }
    
  virtual double giveTransformationJacobian(Domain *d, IntArray &nodes, const FloatArray &lcoords, const IntArray* span, double time) ;
  virtual double giveTransformationJacobian(Domain *d, IntArray &nodes, const FloatArray &lcoords, double time) ;  
  /**
     Returns indices (zero based) of nonzero basis functions for given knot span 
  */
  int giveKnotBasisFuncMask (const IntArray& knotSpan, IntArray& mask) ;
  /** Returns the number of nonzero basis functions at individual knot span */
  int  giveNumberOfKnotBasisFunctions () ;

  /// Returns class name of the receiver.
  const char *giveClassName() const { return "BSplineInterpolation"; }

protected:
  /**
     Evaluates the nonvanishing basis functions of 1d BSpline (algorithm A2.2 from NURBS book) 
     @param i knot span index (zero based) 
     @param u value at which to evaluate
     @param p degree
     @param U knot vector
     @param N computed p nonvanishing functions (N_{i-p,p}-N_{i,p})
  */
  void basisFuns (FloatArray &N, int i, double u, int p, const double* U) ;
  /*
    Compute nonzero basis functions and their derivatives at u of the B-SPline curve
    
    For information on the algorithm, see A2.3 on p72 of 
    the NURBS book. The result is stored in the ders matrix, where 
    ders is of size (n+1,deg+1) and the derivative 
    C'(u) = ders(1,span-deg+j) where j=0...deg
		C''(u)= ders(2,span-deg+j) where j=0...deg
    
    @param n  the degree of the derivation
    @param u  the parametric value
		@param span knot span index (zero based) 
    //@param  span  the span for the basis functions
    @param deg  the degree of the curve
    @param U  the knot vector on which the Basis functions must be computed.
    @param ders  A matrix containing the derivatives of the curve.
    
    @warning n, u and span must be in a valid range.
  */
  void DersBasisFuns(int n, double u, int span, int deg, double* const U, FloatMatrix& ders);
  /*
    Determines the knot span index (Algorithm A2.1 from the NURBS book)

    Determines the knot span for which there exists non-zero basis 
    functions. The span is the index  k for which the parameter 
    u is valid in the [u_k,u_{k+1}] range.
    
    @param n number of control points // HUHU (n+1)
    @param u  the parametric value
    @param p the degree
    @param U knot vector
    @return the span index at u (zero based)
    @warning u must be in a valid range
  */
  int findSpan(int n, int p, double u, const double* U) const;
  /** Return the range of nonzero basis functions for given knot span i and degree p*/
  void giveNonzeroBasisFunctIntervalOnCurve (int i, int p, int &s, int &e) {s=i-p; e=i;}

}; // enf of BSplineInterpolation class definition



/**
 * IntegrationElement represent nonzero knot span, derived from Integration Rule;
 * 
 */
class IGA_IntegrationElement : public GaussIntegrationRule {
protected:
	IntArray knotSpan; // knot_span(nsd)
public:
	IGA_IntegrationElement (int _n, Element* _e, IntArray& _knotSpan) : 
  GaussIntegrationRule (_n, _e, 0, 0, false), 
    knotSpan (_knotSpan) {}
  const IntArray* giveKnotSpan () {return &this->knotSpan;}
  void setKnotSpan(IntArray& src) {this->knotSpan=src;}
};



/**
   Implements base IGAElement, supposed to be a parent class of all elements with B-spline or NURBS based interpolation.
 */
class IGAElement : public Element {
protected:
  // FEInterpolation interpolation; 
  //BSplineInterpolation* interpolation; // HUHU
public:

  IGAElement(int n, Domain *aDomain) : Element (n, aDomain) {}

  IRResultType initializeFrom(InputRecord *ir);
protected:
  virtual int giveNsd() = 0;
};


        
/**
 * StructuralElementEvaluator - base class of all structural elements
 * Individual elements supposed to be derived from StructuralElementEvaluator and IGAElement
 */
class StructuralElementEvaluator {
 protected:
  FloatMatrix* rotationMatrix;
 public:
  StructuralElementEvaluator ();
  virtual void     giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);

 protected:
  virtual Element* giveElement()=0;
  virtual void computeNMatrixAt (FloatMatrix& answer, GaussPoint* gp)=0;
  virtual void computeBMatrixAt (FloatMatrix& answer, GaussPoint* gp)=0;
  virtual void computeStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode, TimeStep *tStep);
  virtual double computeVolumeAround(GaussPoint *gp) { return 0.; }
  /**
   * Updates rotation matrix r(l)=T r(g*) between  local and global coordinate system
   * taking into account also possible local - coordinate system in some elements
   * nodes.
   * Default implementation uses \ref computeGtoLRotationMatrix and
   \ref computeGNDofRotationMatrix  services to compute result.
   * Default implementation uses cached rotation matrix in
   * rotationMatrix attribute, so rotation matrix is computed only once.
   * @return nonzero if transformation is necessary.
   */
  virtual int   updateRotationMatrix() {return 0;}    // to be moved from structural element
  /**
   * Assembles the code numbers of given integration element (sub-patch)
   * This is done by obtaining list of nonzero shape functions and
   * by collecting the code numbers of nodes corresponding to these 
   * shape functions
   * @returns returns nonzero if integration rule code numbers differ from element code numbers
   */
  virtual int giveIntegrationElementCodeNumbers (IntArray& answer, Element* elem, 
                                                 IntegrationRule* ie, EquationID ut) ;
};


/**
 * general purpose Plane stress structural element evaluator
 */
class PlaneStressStructuralElementEvaluator : public StructuralElementEvaluator {
public:
  PlaneStressStructuralElementEvaluator() : StructuralElementEvaluator() {}

protected:
  /// Cached transformation matrix of receiver
  FloatMatrix *rotationMatrix; // to be moved from structural element


  /** Assemble interpolation matrix at given IP
   *  In case of IGAElements, N is assumed to contain only nonzero interpolation functions
   */
  void computeNMatrixAt (FloatMatrix& answer, GaussPoint* gp);
  /** Assembles the strain-displacement matrix of the receiver at given integration point
   *  In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions
   */ 
  void computeBMatrixAt (FloatMatrix& answer, GaussPoint* gp) ;
  double computeVolumeAround(GaussPoint *gp) ;
}; // end of PlaneStressStructuralElementEvaluator definition



class BsplinePlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator {
protected:
  BSplineInterpolation interpolation;

public:

  BsplinePlaneStressElement(int n, Domain *aDomain);
  IRResultType initializeFrom(InputRecord *ir) {
    IGAElement::initializeFrom(ir);
    //PlaneStressStructuralElementEvaluator::initializeFrom(ir);
    return IRRT_OK;
  }
  void     giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
    PlaneStressStructuralElementEvaluator::giveCharacteristicMatrix (answer, mtrx, tStep);}

  
  FEInterpolation *giveInterpolation() {return &this->interpolation;}
  virtual Element* giveElement() {return this;}
protected:
  virtual int giveNsd() {return 2;}
  

};

#endif //iga_h
