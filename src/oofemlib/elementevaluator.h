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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef elementevaluator_h
#define elementevaluator_h

#include "chartype.h"
#include "valuemodetype.h"

#include "femcmpnn.h"
#include "floatmatrix.h"
#include "floatarray.h"

#include "alist.h"
#include "intarray.h"
#include "error.h"
#include "integrationrule.h"
#include "chartype.h"
#include "elementgeometrytype.h"
#include "equationid.h"
#include "valuemodetype.h"
#include "internalstatetype.h"
#include "internalstatevaluetype.h"
#include "elementextension.h"
#include "entityrenumberingscheme.h"
#include "unknowntype.h"
#include "unknownnumberingscheme.h"

#define _IFT_ElementEvaluator_bodyload "bodyloads"
#define _IFT_ElementEvaluator_boundaryload "boundaryloads"


/**
 * Abstract base class for all element evaluators.
 * This class was created by decouplig original Element class
 * into ElementEvaluator and Element Geometry
 * ElementEvaluator evaluates element matricies, load vectors etc.
 * while ElementGeometry store all the geometrical info about the element
 
 * Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 *
 * In general, element class member functions are called in following cases:
 * - When Engineering model assembles governing equation(s) from element's contributions.
 *   Typically when assembles some characteristic matrix of analyzed structure (method assemble),
 *   asks each element for its code numbers and for corresponding characteristic matrix, vector or
 *   value. This will typically uses above mentioned general method for obtaining characteristic
 *   matrices, vectors and values from elements.
 * - When Engineering model has some characteristic matrix stored in sparse form, then is necessary
 *   to build internal profile of such sparse matrix. Engineering model then typically calls
 *   buildInternalStructure member function of sparse matrix. This method then requests element code
 *   numbers from elements and builds its internal profile. This could look strange, because
 *   class sparse matrix should know "something about elements", but only sparse matrix knows, how
 *   to build its internal structure (this may require one or more loops over elements code number).
 * - When element computes its contribution, then it communicates with its cross-section (and cross
 *   section with corresponding material). For some cross-section or material models some further
 *   communication between these classes and element may be necessary (for example in cases of
 *   layered cross section model strains in each layer has to be evaluated from "integrated"
 *   strains (integrated means element-like strains including curvatures in case of beams and plates)
 *   by element, in cases of non local material and many other examples can be found).
 */

namespace oofem{
class FloatArray;
class FloatMatrix;
class Load;
class TimeStep;
class BoundaryLoad;

class OOFEM_EXPORT ElementEvaluator
{

protected:
	/**
	* Array containing indexes of loads (body loads and boundary loads are kept separately),
	* that apply on receiver.
	*/
	IntArray bodyLoadArray, boundaryLoadArray;

public:
//Constructor    
  ElementEvaluator();
// Destructor.
  ~ElementEvaluator(){};
  


  // Overloaded methods:
  virtual IRResultType initializeFrom(InputRecord *ir);
  virtual void giveInputRecord(DynamicInputRecord &input);
  virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
  virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
  //virtual void printOutputAt(FILE *file, TimeStep *tStep);
  virtual const char *giveClassName() const { return "ElementEvaluator"; }

public:
  /// Returns array containing load numbers of loads acting on element
  IntArray *giveBodyLoadArray();
  /// Returns array containing load numbers of boundary loads acting on element.
  IntArray *giveBoundaryLoadArray();


  /**@name Methods referring to code numbers */
  //@{
  /**
  * Returns the location array (array of code numbers) of receiver for given numbering scheme.
  * Results are cached at receiver for default scheme in locationArray attribute.
  */
  virtual void giveLocationArray(IntArray & locationArray, EquationID, const UnknownNumberingScheme & s, ElementGeometry* elementGeometry, IntArray * dofIds = NULL) const;
  virtual void giveLocationArray(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, ElementGeometry* elementGeometry, IntArray *dofIds = NULL) const;
  /**
  * Returns the location array for the boundary of the element.
  * Only takes into account nodes in the bNodes vector.
  */
  virtual void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, EquationID eid, const UnknownNumberingScheme &s, ElementGeometry* elementGeometry, IntArray *dofIds = NULL);
  virtual void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const IntArray &dofIDMask, const UnknownNumberingScheme &s, ElementGeometry* elementGeometry, IntArray *dofIds = NULL);
  /**
  * @return Number of DOFs in element.
  */
  virtual int giveNumberOfDofs() { return 0; }
 
 

  //@}

  /**
  * Returns transformation matrix from global c.s. to local element
  * c.s., i.e. @f$ r_l =T r_g @f$.
  * If no transformation is necessary then answer is empty matrix and zero value is returned.
  * @param answer Computed rotation matrix.
  * @return Nonzero if transformation is necessary, zero otherwise.
  */
  virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
  /**
  * Transformation matrices updates rotation matrix between element-local and primary DOFs,
  * taking into account nodal c.s. and master DOF weights.
  * @param answer Contains the rotation matrix on exit.
  * @param eid Equation ID to compute rotation matrix for.
  * @return True if there is a rotation required, false otherwise.
  */
  virtual bool giveRotationMatrix(FloatMatrix &answer, EquationID eid, ElementGeometry* elementGeometry);
  /**
  * Returns transformation matrix for DOFs from global coordinate system
  * to local coordinate system in nodes.
  * Also includes transformations to slave DOFs.
  * If no transformation is necessary sets answer to empty matrix and returns false.
  * Local stiffness matrix of element should be rotated with answer before assembly.
  * @note Function does most likely NOT need to be overridden.
  * @param answer Computed rotation matrix.
  * @param nodes Nodes to include in element local ordering.
  * @param includeInternal Determines whether or not to include internal dof managers.
  * @param eid Equation ID.
  * @return True if transformation is necessary, false otherwise.
  */
  virtual bool computeDofTransformationMatrix(FloatMatrix &answer, const IntArray &nodes, bool includeInternal, EquationID eid, ElementGeometry* elementGeometry);



  /**@name General element functions */
  //@{
  /**
  * Returns local vector of unknowns. Local vector of unknowns is extracted from
  * engineering model global unknowns vector (if specific dof has assigned
  * corresponding equation number) and from boundary conditions (if dof has active boundary
  * and possibly initial condition). Because unknowns are obtained from engineering
  * model, this must support queries for given unknown and unknown mode. Also
  * engineering model must not be able to complete request for any TimeStep, because
  * keeping history of all unknowns may became impossible. But generally all engineering models
  * should be able return supported unknowns at current and previous time step. Consult
  * reference manual for particular engineering model.
  *
  * @param type   Identifies unknown type (eg. displacement or temperature vector).
  * @param u      Identifies mode of unknown (eg. total value or velocity of unknown).
  * @param tStep  Time step, when vector of unknowns is requested.
  * @param answer Local vector of unknowns.
  */
  virtual void computeVectorOf(EquationID type, ValueModeType u, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry);
  /**
  * Boundary version of computeVectorOf.
  * @param bNodes Boundary nodes.
  * @param eid Equation ID for unknowns.
  * @param u Identifies mode of unknown (eg. total value or velocity of unknown).
  * @param tStep Time step, when vector of unknowns is requested.
  * @param answer Local vector of unknowns.
  */
  virtual void computeBoundaryVectorOf(const IntArray &bNodes, EquationID eid, ValueModeType u, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry);
  /**
  * Returns local vector of unknowns. Local vector of unknowns is extracted from
  * given field and from boundary conditions (if dof has active boundary
  * and possibly initial condition). Because unknowns are obtained from given field
  * model, this must support queries for given unknown.
  *
  * @param field  Source field (eg. displacement or temperature vector).
  * @param u      Value mode of unknown (incremental, total, ...).
  * @param tStep  Time step, when vector of unknowns is requested.
  * @param answer Local vector of unknowns.
  */
  virtual void computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry);
  /**
  * Returns local vector of prescribed unknowns. Local vector of prescribed unknowns is
  * extracted from nodal (and side - if they hold unknowns) boundary conditions.
  *
  * @param ut     Identifies mode of unknown (eg. total values or velocity of unknown).
  * @param type   Value mode of unknown (incremental, total, ...).
  * @param tStep  Time step, when vector of prescribed unknowns is requested.
  * @param answer Local vector of prescribed unknowns. If unknown is not prescribed,
  * zero value is placed on position of free dof.
  */
  void computeVectorOfPrescribed(EquationID ut, ValueModeType type, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry);

  /**
  * Computes or simply returns total number of element's local DOFs.
  * Must be defined by particular element.
  * @return Number of local DOFs of element.
  */
  virtual int computeNumberOfDofs() { return 0; }
  /**
  * Computes the total number of element's global dofs.
  * The transitions from global c.s. to nodal c.s. should NOT be included.
  * @return Total number of global DOFs of element.
  */
  virtual int computeNumberOfGlobalDofs();
  /**
  * Computes the total number of element's primary master DOFs.
  * @param eid ID of equation that DOFs belong to.
  * @return Total number of DOFs belonging to eid.
  */
  int computeNumberOfPrimaryMasterDofs(EquationID eid, ElementGeometry *elementGeometry);
  // initialization to state given by initial conditions
  /** Initialization according to state given by initial conditions.
  * Some type of problems may require initialization of state variables
  * stored in integration points (in statuses related to material models)
  * according to initial conditions. Default implementation is not provided.
  * Typically, loop over all integration points, and call to some initialization
  * method of material (with necessary arguments) have to be made.
  * @param timeStepWhenICApply Time step when IC applies.
  */
  virtual void initializeYourself(TimeStep *timeStepWhenICApply) { }


  /**
  * Tests if the element implements required extension. ElementExtension type defines
  * the list of all available element extensions.
  * @param ext Extension to be tested.
  * @return Nonzero if extension supported.
  * @see ElementExtension
  */
  virtual int testElementExtension(ElementExtension ext) { return 0; }
  /**
  * Updates internal element state (in all integration points of receiver)
  * before nonlocal averaging takes place. Used by so nonlocal materials,
  * because their response in particular point depends not only on state in this point, but
  * depends also on state in point's neighborhood. Nonlocal quantity is computed as nonlocal
  * average of local quantities. Therefore, before updating integration point state depending on
  * nonlocal quantity (or quantities), local quantities in all integration points must be updated
  * in advance. This function should be  overloaded by derived analysis-specific
  * class, which updates state in all receiver's integration points (using updateBeforeNonlocalAverage
  * member function declared at corresponding analysis specific material base class)
  * depending on driving variable (for example - strain vector in case of structural-analysis elements).
  * @param tStep Time step.
  */
  virtual void updateBeforeNonlocalAverage(TimeStep *tStep) { }
  

     /**
     * @name General methods for obtaining element contributions
     * Note: These member functions have to  be overloaded by derived analysis-specific
     * evaluators in order to invoke proper method according to type of component requested.
     * Derived classes from these analysis-specific classes should not modify these functions.
     */
    //@{
    /**
     * Computes characteristic matrix of receiver of requested type in given time step.
     * @param answer Requested characteristic matrix (stiffness, tangent, ...).
     * If element has no capability to compute requested type of characteristic matrix
     * error function is invoked.
     * @param type   Id of characteristic component requested.
     * @param tStep  Time step when answer is computed.
     */
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) = 0;
    /**
     * Computes characteristic vector of receiver of requested type in given time step.
     * If element has no capability to compute requested type of characteristic vector
     * error function is invoked.
     * @param answer Requested characteristic vector.
     * @param type   Id  of characteristic component requested.
     * @param mode   Determines mode of answer.
     * @param tStep  Time step when answer is computed.
     */
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) = 0;
    /**
     * Computes characteristic value of receiver of requested type in given time step.
     * If element has no capability to compute requested type of characteristic value
     * error function is invoked.
     * @param type  Id of characteristic component requested.
     * @param tStep Time step when answer is computed.
     * @return Requested value.
     */
    virtual double giveCharacteristicValue(CharType type, TimeStep *tStep) = 0;
    //@}

    /**
     * @name General methods for computing the contribution from loads
     */
    //@{
    /**
     * Computes the contribution of the given load.
     * @param answer Requested contribution of load.
     * @param load   Load to compute contribution from.
     * @param type   Type of the contribution.
     * @param mode   Determines mode of answer.
     * @param tStep  Time step when answer is computed.
     */
    virtual void computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep) = 0;
    /**
     * Computes the contribution of the given load at the given boundary.
     * @note Elements which do not have an contribution should resize the vector to be empty.
     * @param answer Requested contribution of load.
     * @param load Load to compute contribution from.
     * @param boundary Boundary number.
     * @param type Type of the contribution.
     * @param mode Determines mode of answer.
     * @param tStep Time step when answer is computed.
     */
    virtual void computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep) = 0;
    /**
     * Computes the contribution of the given load at the given boundary edge.
     * @note Elements which do not have an contribution should resize the vector to be empty.
     * @param answer Requested contribution of load.
     * @param load Load to compute contribution from.
     * @param edge Edge number.
     * @param type Type of the contribution.
     * @param mode Determines mode of answer.
     * @param tStep Time step when answer is computed.
     */
    virtual void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep) = 0;
    //@}




};


class EvaluatorMatrix
{
 protected:
  ElementEvaluator** values;
  /// Number of rows.
  int nRows;
  /// Number of columns.
  int nColumns;
  
 public:
  
 EvaluatorMatrix(int n, int m)
    {
      
      values = new ElementEvaluator*[n*m];
    }
  
 EvaluatorMatrix(int n)
    {
      
      values = new ElementEvaluator*[n*n];
    }
  
 EvaluatorMatrix(int n, int m, ElementEvaluator** eA)
    {
      values = new ElementEvaluator*[n*m];
      for(int i = 0; i< m*n;i++)
	{
	  values[i] = eA[i];
	}
      
    }
  
 EvaluatorMatrix(int n, ElementEvaluator** eA)
    {
      values = new ElementEvaluator*[n*n];
      for(int i = 0; i< n*n;i++)
	{
	  values[i] = eA[i];
	}
      
      
    }
  
 EvaluatorMatrix()
    {
      values = NULL;
    }
  
  
  ~EvaluatorMatrix()
    {
      if(values != NULL)
	{
	  delete [] values;
	  values = NULL;
	}
      
    }
   int giveNumberOfRows() const { return nRows; }
    /// Returns number of columns of receiver.
    int giveNumberOfColumns() const { return nColumns; }
  
  
  ElementEvaluator* at(int i, int j) const {return values [ ( j - 1 ) * nRows + i - 1 ];}
  ElementEvaluator* &at(int i, int j) { return values [ ( j - 1 ) * nRows + i - 1 ]; }
  
  
};


} // end namespace oofem
#endif // elementevaluator_h
