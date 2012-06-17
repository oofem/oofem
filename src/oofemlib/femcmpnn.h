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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef femcmpnn_h
#define femcmpnn_h

#ifndef __MAKEDEPEND
 #include <string>
#endif

#include "error.h"
#include "interfacetype.h"
#include "inputrecord.h"
#include "classtype.h"
#include "entityrenumberingscheme.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

class DataStream;
class Domain;
class Interface;
class TimeStep;

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
class FEMComponent
{
protected:
    /// Component number
    int number;
    /// Link to domain object, useful for communicating with other FEM components
    Domain *domain;

public:
    /// Constructor used for instanciating temporary objects.
    FEMComponent() { }
    /**
     * Regular constructor, creates component with given number and belonging to given domain.
     * @param n Component number in particular domain. For instance, can represent
     * node number in particular domain.
     * @param d Domain to which component belongs to.
     */
    FEMComponent(int n, Domain *d) {
        number = n;
        domain = d;
    }
    /// Virtual destructor.
    virtual ~FEMComponent() { }

    /**
     * Returns classType id of receiver. Intended for run time
     * type checking. Every derived class have to overload this method.
     * @see classType.
     * @return Class type of receiver.
     */
    virtual classType giveClassID() const { return FEMComponentClass; }
    /// @return Class name of the receiver.
    virtual const char *giveClassName() const  = 0;
    /// @return Input record name of the receiver.
    virtual const char *giveInputRecordName() const { return ( this->giveClassName() ); }
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
    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    /**
     * Setups the input record string of receiver.
     * @param str String to be filled by input record.
     * @param keyword Determines if record keyword should be printed.
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Stores receiver state to output stream.
     * Writes the FEMComponent class-id in order to allow test whether correct data are then restored.
     * @param stream Output stream.
     * @param mode Determines amount of info required in stream (state, definition, ...).
     * @param obj Special parameter, used only to send particular integration point to material class version of this method.
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * Reads the FEMComponent class-id in order to allow test consistency.
     * @see saveContext
     * @param stream Input stream.
     * @param mode Determines amount of info available in stream (state, definition, ...).
     * @param obj Special parameter for sending extra information.
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Allows programmer to test some internal data, before computation begins.
     * For example, one may use this function, to ensure that element has material with
     * required capabilities is assigned to element. This must be done after all
     * mesh components are instanciated.
     * @return Nonzero if receiver is consistent.
     */
    virtual int checkConsistency() { return 1; }
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

    /**
     * @name error and warning reporting methods
     * These methods will print error (or warning) message using oofem default loggers.
     * Do not use these methods directly, to avoid specify file and line parameters.
     * More preferably, use these methods via corresponding OOFEM_CLASS_ERROR and OOFEM_CLASS_WARNING macros,
     * that will include file and line parameters automatically.
     *
     * Uses variable number of arguments, so a format string followed by optional arguments is expected
     * (according to printf conventions).
     * @param file  source file name, where error encountered (where error* function called)
     * @param line  source file line number, where error encountered
     */
    //@{
    /// Prints error message and exits.
    void error(const char *file, int line, const char *format, ...) const;
    /// Prints warning message.
    void warning(const char *file, int line, const char *format, ...) const;
    //@}
#ifdef __OOFEG
    virtual void drawYourself(oofegGraphicContext &) { }
#endif
};
} // end namespace oofem
#endif // femcmpnn_h








