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
  oofem nodes - control points (coordinates ) + dofs
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


#ifdef __OOFEG
class StructuralElementEvaluator;
void drawIGAPatchDeformedGeometry(Element* elem, StructuralElementEvaluator* se, oofegGraphicContext &gc, UnknownType);
#endif




class FEIIGAElementGeometryWrapper : public FEICellGeometry 
{
public:
  const IntArray *knotSpan;
  Element* elem;
public:
  FEIIGAElementGeometryWrapper (Element* _elem, const IntArray* _knotSpan) : FEICellGeometry() {
    this->elem = _elem; this->knotSpan = _knotSpan; }
  FEIIGAElementGeometryWrapper (Element* _elem) : FEICellGeometry() {
    this->elem = _elem; this->knotSpan = NULL; }

  int giveNumberOfVertices () const {return elem->giveNumberOfNodes();}
  const FloatArray* giveVertexCoordinates(int i) const {return elem->giveNode(i)->giveCoordinates();}
};



class BSplineInterpolation : public FEInterpolation  
{
protected:
	/// number of spatial directions
	int nsd;        
	/// degree in each direction 
	int *degree;                               // eg. 2
	/// knotValues[nsd]
	FloatArray *knotValues;                    // eg. 0 1 2 3 4 5
	/// knotMultiplicity[nsd]
	IntArray *knotMultiplicity;                // eg. 3 1 1 1 2 3
	/// numberOfControlPoints[nsd]
	/// for TSpline this is filled by values corresponding to case when there are no T-junctions
	/// (i.e. TSpline = BSpline)
	int *numberOfControlPoints; 
	/// knotVectors[nsd][knot_vector_size]
	double **knotVector;                       // eg. 0 0 0 1 2 3 4 4 5 5 5
	/// nonzero spans in each directions[nsd]
	int *numberOfKnotSpans;                    // eg. 5 (0-1,1-2,2-3,3-4,4-5)
public:
  BSplineInterpolation (int nsd) : FEInterpolation (0) {this->nsd=nsd;}
  ~BSplineInterpolation();
	/**
	 * Returns number of spatial dimensions
	 */
	int const giveNsd() {return nsd;}
  IRResultType initializeFrom(InputRecord *ir);                   
  virtual int giveNumberOfKnotSpans(int dim) {return numberOfKnotSpans[dim-1];}
  virtual int giveNumberOfControlPoints(int dim) {return numberOfControlPoints[dim-1];}
  virtual double** const giveKnotVector() {return this->knotVector;}
  virtual IntArray* const giveKnotMultiplicity(int dim) {return &this->knotMultiplicity[dim-1];}
  virtual FloatArray* const giveKnotValues(int dim) {return &this->knotValues[dim-1];}
  /**
   * Evaluates the array of interpolation functions (shape functions) at given point.
   * @param answer contains resulting array of evaluated interpolation functions
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   *
   * see also giveNonzeroBasisFunctMask method
   */
  virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  /**
   * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
   * These derivatives are in global coordinate system (where the nodal coordinates are defined)
   * @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   */
  virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  /**
   * Evaluates global coordinates from given local ones
   * @param answer contains resulting global coordinates
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   */
  virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
    OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
    return 0;
  }
  /**
   * Evaluates the jacobian of transformation between local and global coordinates.
   */
  virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;


  /**
     Returns indices (zero based) of nonzero basis functions for given knot span 
  */
  int giveKnotSpanBasisFuncMask (const IntArray& knotSpan, IntArray& mask) ;
  /** Returns the number of nonzero basis functions for given knot span */
  int giveNumberOfKnotSpanBasisFunctions (const IntArray& knotSpan) ;

  /// Returns class name of the receiver.
  const char *giveClassName() const { return "BSplineInterpolation"; }
  virtual bool hasSubPatchFormulation() {return true;}
protected:
  /**
     Evaluates the nonvanishing basis functions of 1d BSpline (algorithm A2.2 from NURBS book) 
     @param span knot span index (zero based) 
     @param u value at which to evaluate
     @param p degree
     @param U knot vector
     @param N computed p+1 nonvanishing functions (N_{span-p,p}-N_{span,p})

		 @warning u and span must be in a valid range.
  */
  void basisFuns (FloatArray &N, int span, double u, int p, const double* U) ;
  /**
    Computes nonzero basis functions and their derivatives at u
    
    For information on the algorithm, see A2.3 on p72 of 
    the NURBS book. The result is stored in the ders matrix, where 
    ders is of size (n+1,p+1) and the derivative 
		N(u)  = ders(0,span-p+j) where j=0...p 
    N'(u) = ders(1,span-p+j) where j=0...p
		N''(u)= ders(2,span-p+j) where j=0...p
    
    @param n the degree of the derivation
    @param u the parametric value
		@param span knot span index (zero based) 
    @param p the degree
    @param U knot vector
    @param ders matrix containing the derivatives of the basis functions.
    
    @warning n, u and span must be in a valid range.
  */
  void dersBasisFuns(int n, double u, int span, int p, double* const U, FloatMatrix& ders);
  /**
    Determines the knot span index (Algorithm A2.1 from the NURBS book)

    Determines the knot span for which there exists non-zero basis 
    functions. The span is the index k for which the parameter 
    u is valid in the (u_k,u_{k+1}] range.
    
    @param n number of control points - 1 (number of ctrl pnts = n + 1)
    @param u the parametric value
    @param p the degree
    @param U knot vector
    @return the span index at u (zero based)
    @warning u must be in a valid range
  */
  int findSpan(int n, int p, double u, const double* U) const;
  /** Return the range of nonzero basis functions for given knot span and given degree */
  void giveNonzeroBasisFuncInterval (int span, int deg, int &s, int &e) {s=span-deg; e=span;}

}; // enf of BSplineInterpolation class definition



class NURBSInterpolation : public BSplineInterpolation  
{
public:
  NURBSInterpolation (int nsd) : BSplineInterpolation (nsd) {}
  ~NURBSInterpolation();

  /**
   * Evaluates the array of interpolation functions (shape functions) at given point.
   * @param answer contains resulting array of evaluated interpolation functions
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   *
   * see also giveNonzeroBasisFunctMask method of BSplineInterpolation
   */
  virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  /**
   * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
   * These derivatives are in global coordinate system (where the nodal coordinates are defined)
   * @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   */
  virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  /**
   * Evaluates global coordinates from given local ones
   * @param answer contains resulting global coordinates
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   */
  virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
    OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
    return 0;
  }
  /**
   * Evaluates the jacobian of transformation between local and global coordinates.
   */
  virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;

  /// Returns class name of the receiver.
  const char *giveClassName() const { return "NURBSInterpolation"; }
}; // end of NURBSInterpolation class definition



/* alternatively, it is possible to store for individual control points open local knot vector;
	 however, this is not enough as I need to know how many knots have been prepended and appended
	 (see createLocalKnotVector) as these are relevant for finding relevant knot span using BSpline method
	 (this could be overcome by writing corresponding TSpline method) and for extraction proper
	 basis function and its derivatives (computed by BSpline methods) */

class TSplineInterpolation : public BSplineInterpolation  
{
 protected:
	/// localIndexKnotVector[number_of_control_points][nsd][degree+2]
	int ***localIndexKnotVector;
	int totalNumberOfControlPoints;
	/// temporary open local knot vector to enable use of BSpline algorithms (common for all directions)
	/// openLocalKnotVector[3*max_degree+2]
	double *openLocalKnotVector;
public:
  TSplineInterpolation (int nsd) : BSplineInterpolation (nsd) {}
  ~TSplineInterpolation();

  IRResultType initializeFrom(InputRecord *ir);                   
	void setNumberOfControlPoints(int num) {this->totalNumberOfControlPoints=num;}
  /**
   * Evaluates the array of interpolation functions (shape functions) at given point.
   * @param answer contains resulting array of evaluated interpolation functions
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   *
   * see also giveNonzeroBasisFunctMask method of TSplineInterpolation
   */
  virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  /**
   * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
   * These derivatives are in global coordinate system (where the nodal coordinates are defined)
   * @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   */
  virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  /**
   * Evaluates global coordinates from given local ones
   * @param answer contains resulting global coordinates
   * @param lcoords array containing (local) coordinates
   * @param cellgeo underlying cell geometry
   * @param time time
   */
  virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;
  virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) {
    OOFEM_ERROR ("Not yet inplemented, contact lazy dr for implementation");
    return 0;
  }
  /**
   * Evaluates the jacobian of transformation between local and global coordinates.
   */
  virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry& cellgeo, double time) ;

  /**
     Returns indices (zero based) of nonzero basis functions for given knot span 
  */
  int giveKnotSpanBasisFuncMask (const IntArray& knotSpan, IntArray& mask) ;
  /** Returns the number of nonzero basis functions for given knot span */
  int giveNumberOfKnotSpanBasisFunctions (const IntArray& knotSpan) ;

  /// Returns class name of the receiver.
  const char *giveClassName() const { return "TSplineInterpolation"; }
protected:
  /**
     Evaluates the middle basis function on local knot vector at u
     @param u value at which to evaluate
     @param p degree
     @param U global knot values
		 @param I local index knot vector
     @return N computed middle basis function

		 @warning u must be in a valid range.
  */
  double basisFunction (double u, int p, const FloatArray& U, const int* I) ;
  /*
    Computes the middle basis function and it derivatives on local knot vector at u

    The result is stored in the ders vector
		N(u)  = ders(0)
    N'(u) = ders(1)
		N''(u)= ders(2)
    
    @param n the degree of the derivation
    @param u the parametric value
    @param p the degree
		@param U global knot values
		@param I local index knot vector
    @param ders vector containing the derivatives of the middle basis function.
    
    @warning n and u must be in a valid range.
  */
  void dersBasisFunction(int n, double u, int p, const FloatArray& U, const int* I, FloatArray& ders);
  /*
    Creates local open knot vector.
		This is generally done extracting knot values from global knot vector using the local index knot vector
		and by prepending p times the first knot and appending p times the last knot.
		However, existing knot multiplicity at the start and end must be accounted for.
		
    @param p the degree
    @param U global knot values
    @param I local index knot vector
    @param prepent number of prepended entries
		@param append number of appended entries
  */
	void createLocalKnotVector(int p, const FloatArray& U, const int* I, int* prepend, int* append);
  /**
     Returns indices (zero based) of nonzero basis functions for given knot span interval (from start to end)
  */
  int giveKnotSpanBasisFuncMask (const IntArray& startKnotSpan, const IntArray& endKnotSpan, IntArray& mask) ;
  /** Returns the number of nonzero basis functions at given knot span interval (from start to end) */
  int  giveNumberOfKnotSpanBasisFunctions (const IntArray& startKnotSpan, const IntArray& endKnotSpan) ;
}; // end of TSplineInterpolation class definition


/**
 * IntegrationElement represent nonzero knot span, derived from Integration Rule;
 * 
 */
class IGAIntegrationElement : public GaussIntegrationRule {
protected:
	IntArray knotSpan; // knot_span(nsd)
public:
	IGAIntegrationElement (int _n, Element* _e, IntArray& _knotSpan) : 
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
public:
  IGAElement(int n, Domain *aDomain) : Element (n, aDomain) {}
  IRResultType initializeFrom(InputRecord *ir);

#ifdef __OOFEG
    //
    // Graphics output
    //
  virtual void  drawRawGeometry(oofegGraphicContext &mode);
#endif


protected:
  virtual int giveNsd() = 0;   // this info is available also from interpolation. Do we need it here ???
};


/**
	 IGATSplineElement setups integration rules differently from IGAElement
*/
class IGATSplineElement : public IGAElement {
public:
  IGATSplineElement(int n, Domain *aDomain) : IGAElement (n, aDomain) {}
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
  virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
  virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) {
    if ( type == ElementNonForceLoadVector )  {
        this->computeNonForceLoadVector(answer, tStep, mode);
    } else {
      answer.resize(0);
    }
  }

 protected:
  virtual Element* giveElement()=0;
  virtual void computeNMatrixAt (FloatMatrix& answer, GaussPoint* gp)=0;
  virtual void computeBMatrixAt (FloatMatrix& answer, GaussPoint* gp)=0;
  virtual void computeStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode, TimeStep *tStep);
  virtual double computeVolumeAround(GaussPoint *gp) { return 0.; }
  void  computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
  void  computeBcLoadVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
  virtual void giveInternalForcesVector(FloatArray &answer,
                                        TimeStep *, int useUpdatedGpRecord = 0) {
    answer.resize(0);
  }
  void computeVectorOf(EquationID type, ValueModeType u,
											 TimeStep *stepN, FloatArray &answer) {
    this->giveElement()->computeVectorOf(type, u, stepN, answer);}
  void computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *stepN, FloatArray &answer) {
    this->giveElement()->computeVectorOf(field, u, stepN, answer);}
  void computeVectorOfPrescribed(EquationID ut, ValueModeType type, TimeStep *stepN, FloatArray &answer) {
    this->giveElement()->computeVectorOfPrescribed (ut, type, stepN, answer);
  }
  bool   isActivated(TimeStep *atTime) {return true;}
  void   updateInternalState(TimeStep *stepN);
  void   computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
  void   computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
  /* Optimized version, allowing to pass element displacents as parameter.
     Standart version has a huge performance leak; in typical iga element the element vector is VERY large
     and its querying for each point take more time than strain evaluation. And this has to be done for each 
     integration point. This optimized version allows to assemble displacement vector only once (for all IP) 
     and pass this vector as parameter
 */
  void   computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray& u);
  void   computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray& u);

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
  virtual int updateRotationMatrix() {return 0;}    // to be moved from structural element
  /**
   * Assembles the code numbers of given integration element (sub-patch)
   * This is done by obtaining list of nonzero shape functions and
   * by collecting the code numbers of nodes corresponding to these 
   * shape functions
   * @returns returns nonzero if integration rule code numbers differ from element code numbers
   */
  virtual int giveIntegrationElementCodeNumbers (IntArray& answer, Element* elem, 
                                                 IntegrationRule* ie, EquationID ut) ;

  /**
   * Assembles the local element code numbers of given integration element (sub-patch)
   * This is done by obtaining list of nonzero shape functions and
   * by collecting the code numbers of nodes corresponding to these 
   * shape functions
   * @returns returns nonzero if integration rule code numbers differ from element code numbers
   */
  virtual int giveIntegrationElementLocalCodeNumbers (IntArray& answer, Element* elem, 
                                                      IntegrationRule* ie, EquationID ut) ;
#ifdef __OOFEG
  friend void drawIGAPatchDeformedGeometry(Element* elem, StructuralElementEvaluator* se, oofegGraphicContext &gc, UnknownType);
#endif
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
  void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
    answer.resize(2);
    answer.at(1) = D_u;
    answer.at(2) = D_v;
  }

}; // end of PlaneStressStructuralElementEvaluator definition


/**
 * general purpose 3d structural element evaluator
 */
class Space3dStructuralElementEvaluator : public StructuralElementEvaluator {
public:
  Space3dStructuralElementEvaluator() : StructuralElementEvaluator() {}

protected:

  /** Assemble interpolation matrix at given IP
   *  In case of IGAElements, N is assumed to contain only nonzero interpolation functions
   */
  void computeNMatrixAt (FloatMatrix& answer, GaussPoint* gp);
  /** Assembles the strain-displacement matrix of the receiver at given integration point
   *  In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions
   */ 
  void computeBMatrixAt (FloatMatrix& answer, GaussPoint* gp) ;
  double computeVolumeAround(GaussPoint *gp) ;
  void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
    answer.resize(3);
    answer.at(1) = D_u;
    answer.at(2) = D_v;
    answer.at(3) = D_w;
  }

}; // end of SpaceStructuralElementEvaluator definition


class BsplinePlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator {
protected:
  BSplineInterpolation interpolation;

public:
  BsplinePlaneStressElement(int n, Domain *aDomain);
  IRResultType initializeFrom(InputRecord *ir);

  void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
    PlaneStressStructuralElementEvaluator::giveCharacteristicMatrix (answer, mtrx, tStep);}
  virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
    PlaneStressStructuralElementEvaluator::giveCharacteristicVector(answer, type, mode, t);}
  
  FEInterpolation *giveInterpolation() {return &this->interpolation;}
  virtual Element* giveElement() {return this;}
  void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
    PlaneStressStructuralElementEvaluator::giveDofManDofIDMask(inode, u, answer);
  }
  virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans*2; }
  void updateInternalState(TimeStep *stepN) {PlaneStressStructuralElementEvaluator::updateInternalState (stepN);}
#ifdef __OOFEG
    //
    // Graphics output
    //
  virtual void  drawScalar(oofegGraphicContext &context);
	virtual void drawDeformedGeometry(oofegGraphicContext &mode, UnknownType ut) {
		drawIGAPatchDeformedGeometry(this, this, mode, ut);
	}
#endif

protected:
  virtual int giveNsd() {return 2;}
};


class NURBSPlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator {
protected:
  NURBSInterpolation interpolation;

public:
  NURBSPlaneStressElement(int n, Domain *aDomain);
  IRResultType initializeFrom(InputRecord *ir);

  void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
    PlaneStressStructuralElementEvaluator::giveCharacteristicMatrix (answer, mtrx, tStep);}
  virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
    PlaneStressStructuralElementEvaluator::giveCharacteristicVector(answer, type, mode, t);}
  
  FEInterpolation *giveInterpolation() {return &this->interpolation;}
  virtual Element* giveElement() {return this;}
  void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
    PlaneStressStructuralElementEvaluator::giveDofManDofIDMask(inode, u, answer);
  }
  virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans*2; }
  void updateInternalState(TimeStep *stepN) {PlaneStressStructuralElementEvaluator::updateInternalState (stepN);}
#ifdef __OOFEG
    //
    // Graphics output
    //
  virtual void  drawScalar(oofegGraphicContext &context);
	virtual void drawDeformedGeometry(oofegGraphicContext &mode, UnknownType ut) {
		drawIGAPatchDeformedGeometry(this, this, mode, ut);
	}

#endif

protected:
  virtual int giveNsd() {return 2;}
};



class TSplinePlaneStressElement : public IGATSplineElement, public PlaneStressStructuralElementEvaluator {
protected:
  TSplineInterpolation interpolation;

public:
  TSplinePlaneStressElement(int n, Domain *aDomain);
  IRResultType initializeFrom(InputRecord *ir) {
		IGATSplineElement::initializeFrom(ir);
		//PlaneStressStructuralElementEvaluator::initializeFrom(ir);
		return IRRT_OK;
	}

  void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
    PlaneStressStructuralElementEvaluator::giveCharacteristicMatrix (answer, mtrx, tStep);}
  virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
    PlaneStressStructuralElementEvaluator::giveCharacteristicVector(answer, type, mode, t);}
  
  FEInterpolation *giveInterpolation() {return &this->interpolation;}
  virtual Element* giveElement() {return this;}
  void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
    PlaneStressStructuralElementEvaluator::giveDofManDofIDMask(inode, u, answer);
  }
  virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans*2; }
  void updateInternalState(TimeStep *stepN) {PlaneStressStructuralElementEvaluator::updateInternalState (stepN);}
#ifdef __OOFEG
    //
    // Graphics output
    //
  virtual void  drawScalar(oofegGraphicContext &context);
#endif

protected:
  virtual int giveNsd() {return 2;}
};



class NURBSSpace3dElement : public IGAElement, public Space3dStructuralElementEvaluator {
protected:
  NURBSInterpolation interpolation;

public:
  NURBSSpace3dElement(int n, Domain *aDomain);
  IRResultType initializeFrom(InputRecord *ir);

  void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
    Space3dStructuralElementEvaluator::giveCharacteristicMatrix (answer, mtrx, tStep);}
  virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
    Space3dStructuralElementEvaluator::giveCharacteristicVector(answer, type, mode, t);}
  
  FEInterpolation *giveInterpolation() {return &this->interpolation;}
  virtual Element* giveElement() {return this;}
  void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
    Space3dStructuralElementEvaluator::giveDofManDofIDMask(inode, u, answer);
  }
  virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans*3; }
  void updateInternalState(TimeStep *stepN) {Space3dStructuralElementEvaluator::updateInternalState (stepN);}
#ifdef __OOFEG
    //
    // Graphics output
    //
  virtual void  drawScalar(oofegGraphicContext &context) ;
	virtual void drawDeformedGeometry(oofegGraphicContext &mode, UnknownType ut) {
		drawIGAPatchDeformedGeometry(this, this, mode, ut);
	}
#endif

protected:
  virtual int giveNsd() {return 3;}
};


#endif //iga_h
