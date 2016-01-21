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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef femcmpnn_h
#define femcmpnn_h

#include <string>

#include "oofemcfg.h"
#include "interfacetype.h"
#include "irresulttype.h"
#include "contextioresulttype.h"
#include "contextmode.h"

namespace oofem {
class DataStream;
class Domain;
class Interface;
class TimeStep;
class InputRecord;
class DynamicInputRecord;
class oofegGraphicContext;
class EntityRenumberingFunctor;
class FloatArray;
class IntArray;
class FloatMatrix;

/**
 * The top abstract class of all classes constituting the finite element mesh.
 * Defines the attributes and methods common to all components of mesh:
 * elements, nodes, time steps, materials, loads and load-time functions.
 * This class defines the two attributes common to all component classes ;
 * 'number' is primarily used for reading data in the data file. 'domain' is
 * used for communicating with other components (e.g., for an element to obtain its material),
 * for accessing the linear system and the data file.
 * @see error handles error reporting.
 * @see checkConsistency to ensure, whether internal data structures are consistent.
 */
class OOFEM_EXPORT FEMComponent
{
protected:
    /// Component number
    int number;
    /// Link to domain object, useful for communicating with other FEM components
    Domain *domain;

public:
    /**
     * Regular constructor, creates component with given number and belonging to given domain.
     * @param n Component number in particular domain. For instance, can represent
     * node number in particular domain.
     * @param d Domain to which component belongs to.
     */
    FEMComponent(int n, Domain * d) : number(n), domain(d) { }
    /// Virtual destructor.
    virtual ~FEMComponent() { }

    /// @return Class name of the receiver.
    virtual const char *giveClassName() const = 0;
    /// @return Input record name of the receiver.
    virtual const char *giveInputRecordName() const = 0;
    /// @return Domain which receiver belongs to.
    Domain *giveDomain() const { return domain; }
    /**
     * Sets associated Domain
     * @param d New domain which receiver should belong to.
     */
    virtual void setDomain(Domain *d) { this->domain = d; }
    /// @return Component number of receiver.
    int giveNumber() const { return number; }
    /**
     * Sets number of receiver.
     * @param num New number of receiver.
     */
    void setNumber(int num) { this->number = num; }
    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various FEM components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflect updated numbering. The renumbering function is passed, which is supposed
     * to return an updated number of specified entity type based on old number.
     */
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f) { }
    /**
     * Initializes receiver according to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. Input record can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see IR_GIVE_FIELD
     * @see IR_GIVE_OPTIONAL_FIELD
     * @param ir Input record to initialize from.
     * @return IRResultType
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Setups the input record string of receiver.
     * @param input Dynamic input record to be filled by receiver.
     */
    virtual void giveInputRecord(DynamicInputRecord &input);
    /**
     * Stores receiver state to output stream.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @param obj Special parameter, used only to send particular integration point to material class version of this method.
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * @see saveContext
     * @param stream Input stream.
     * @param mode Determines amount of info available in stream (state, definition, ...).
     * @param obj Special parameter for sending extra information.
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /**
     * Allows programmer to test some internal data, before computation begins.
     * For example, one may use this function, to ensure that element has material with
     * required capabilities is assigned to element. This must be done after all
     * mesh components are instanciated.
     * @return Nonzero if receiver is consistent.
     */
    virtual int checkConsistency();
    /**
     * Prints output of receiver to stream, for given time step.
     * This is used for output into the standard output file.
     * @param file File pointer to print to.
     * @param tStep Time step to write for.
     */
    virtual void printOutputAt(FILE *file, TimeStep *tStep) { }
    /// Prints receiver state on stdout. Useful for debugging.
    virtual void printYourself() { }
    /**
     * Interface requesting service.
     * @see InterfaceType
     * @return Requested interface if implemented, otherwise NULL.
     */
    virtual Interface *giveInterface(InterfaceType t) { return NULL; }

    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const;
};
} // end namespace oofem
#endif // femcmpnn_h
