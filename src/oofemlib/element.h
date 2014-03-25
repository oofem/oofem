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

#ifndef element_h
#define element_h

#include "elementgeometry.h"
#include "elementevaluator.h"
#include "baseelement.h"
#include "contextioerr.h"

namespace oofem {



/**
 * Abstract base class for all finite elements. Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * classes in order to invoke proper method according to type of component requested).
 */
	class OOFEM_EXPORT Element : public ElementEvaluator, public ElementGeometry, public BaseElement
{
protected:

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
	Element(int n, Domain* d);
    /// Virtual destructor.
	virtual ~Element(){;}

	virtual ElementGeometry* giveElementGeometry() { return this;}
	virtual ElementEvaluator* giveElementEvaluator() { return this;}
	virtual void printOutputAt(FILE *file, TimeStep *tStep)
	{
		ElementGeometry::printOutputAt(file, tStep);
//		ElementEvaluator::printOutputAt(file, tStep);
	}
	
	IRResultType initializeFrom(InputRecord *ir)
	{
		IRResultType iTeG,iTeE;
		iTeE =  ElementEvaluator::initializeFrom(ir);
		iTeG = ElementGeometry::initializeFrom(ir);
		return iTeE;
	}

	void giveInputRecord(DynamicInputRecord &input)
	{
		ElementEvaluator::giveInputRecord(input);		
		ElementGeometry::giveInputRecord(input);
	}

	contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj)
	{
			contextIOResultType iores;
			if ((iores = ElementGeometry::saveContext(stream, mode, obj)) != CIO_OK) {
				THROW_CIOERR(iores);
			}

			if ((iores = ElementEvaluator::saveContext(stream, mode, obj)) != CIO_OK) {
				THROW_CIOERR(iores);
			}
			return CIO_OK;
	}

	contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj)
	{
		contextIOResultType iores;
		if ((iores = ElementGeometry::restoreContext(stream, mode, obj)) != CIO_OK) {
			THROW_CIOERR(iores);
		}

		if ((iores = ElementEvaluator::restoreContext(stream, mode, obj)) != CIO_OK) {
			THROW_CIOERR(iores);
		}
		return CIO_OK;
	}

	virtual const char *giveClassName() const { return "Element"; }




	virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) {;}
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) {;}
    virtual double giveCharacteristicValue(CharType type, TimeStep *tStep) {return 0;}
    virtual void computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep) {;}
    virtual void computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep) {;}
    virtual void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep) {;}

	/**
	* Returns the location array (array of code numbers) of receiver for given numbering scheme.
	* Results are cached at receiver for default scheme in locationArray attribute.
	*/
	virtual void giveLocationArray(IntArray & locationArray, EquationID id, const UnknownNumberingScheme & s, IntArray * dofIds = NULL)
	{
		ElementEvaluator :: giveLocationArray(locationArray, id, s, this, dofIds);
	}
	virtual void giveLocationArray(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIds = NULL)
	{
		ElementEvaluator::giveLocationArray(locationArray, dofIDMask, s, this, dofIds);
	}
	/**
	* Returns the location array for the boundary of the element.
	* Only takes into account nodes in the bNodes vector.
	*/
	virtual void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, EquationID eid, const UnknownNumberingScheme &s, IntArray *dofIds = NULL)
	{
		ElementEvaluator :: giveBoundaryLocationArray(locationArray, bNodes, eid, s, this, dofIds);
	}
	virtual void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIds = NULL)
	{
		ElementEvaluator::giveBoundaryLocationArray(locationArray, bNodes, dofIDMask, s, this, dofIds);
	}


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
	virtual void computeVectorOf(EquationID type, ValueModeType u, TimeStep *tStep, FloatArray &answer)
	{
		ElementEvaluator :: computeVectorOf(type, u, tStep, answer,this);
	}
	
	/**
	* Boundary version of computeVectorOf.
	* @param bNodes Boundary nodes.
	* @param eid Equation ID for unknowns.
	* @param u Identifies mode of unknown (eg. total value or velocity of unknown).
	* @param tStep Time step, when vector of unknowns is requested.
	* @param answer Local vector of unknowns.
	*/
	virtual void computeBoundaryVectorOf(const IntArray &bNodes, EquationID eid, ValueModeType u, TimeStep *tStep, FloatArray &answer)
	{
		ElementEvaluator :: computeBoundaryVectorOf(bNodes, eid, u, tStep, answer, this);
	}
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
	virtual void computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *tStep, FloatArray &answer)
	{
		ElementEvaluator :: computeVectorOf(field, u, tStep, answer, this);

	}
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
	void computeVectorOfPrescribed(EquationID ut, ValueModeType type, TimeStep *tStep, FloatArray &answer)
	{
		ElementEvaluator :: computeVectorOfPrescribed(ut, type, tStep, answer, this);
	}

	/**
	* Returns transformation matrix from global c.s. to local element
	* c.s., i.e. @f$ r_l =T r_g @f$.
	* If no transformation is necessary then answer is empty matrix and zero value is returned.
	* @param answer Computed rotation matrix.
	* @return Nonzero if transformation is necessary, zero otherwise.
	*/
	virtual bool computeGtoLRotationMatrix(FloatMatrix &answer)
	{
		return ElementEvaluator :: computeGtoLRotationMatrix(answer, this);
	}
	/**
	* Transformation matrices updates rotation matrix between element-local and primary DOFs,
	* taking into account nodal c.s. and master DOF weights.
	* @param answer Contains the rotation matrix on exit.
	* @param eid Equation ID to compute rotation matrix for.
	* @return True if there is a rotation required, false otherwise.
	*/
	virtual bool giveRotationMatrix(FloatMatrix &answer, EquationID eid)
	{
		return ElementEvaluator :: giveRotationMatrix(answer, eid, this);
	}
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
	virtual bool computeDofTransformationMatrix(FloatMatrix &answer, const IntArray &nodes, bool includeInternal, EquationID eid)
	{
		return ElementEvaluator :: computeDofTransformationMatrix(answer, nodes, includeInternal, eid, this);
	}



};
    
} // end namespace oofem
#endif //element_h
