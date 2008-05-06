/* $Header: /home/cvs/bp/oofem/oofemlib/src/femcmpnn.h,v 1.15.4.1 2004/04/05 15:19:43 bp Exp $ */
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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

//
//  CLASS FEMCOMPONENT
//

#ifndef femcmpnn_h
#define femcmpnn_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string>
#endif

#include "intarray.h"
#include "flotarry.h"
#include "domain.h"
#include "error.h"
#include "interface.h"
#include "inputrecord.h"
#include "classtype.h"
#include "interfacetype.h"
#include "entityrenumberingscheme.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

class DataStream;
/**
 * The top abstract class of all classes constituting the finite element mesh.
 * Defines the atributes and methods common to all components of mesh.
 */
class FEMComponent
{
    /*
     * This class is an abstract class, the superclass of all classes that imple-
     * ment the components of a finite element mesh : elements, nodes, time steps,
     * materials, loads and load-time functions.
     * DESCRIPTION
     * This class defines the two attributes common to all component classes ;
     * 'number' is primarily used for reading data in the data file. 'domain' is
     * used for communicating with other components (e.g., for an element to ob-
     * tain its material), for accessing the linear system and the data file.
     * Member function error handles error reporting.
     * Member function checkConsistency() to ensure, whether internal
     * data structures are consistent.
     * TASKS
     *
     */
protected:
    /// component number
    int number;
    /** link to domain object, usefull for communicating
     *  with other FEM components */
    Domain *domain;

public:
    /// Constructor used for instanciating temporary objects.
    FEMComponent() { }                                      // constructors
    /** Regular constructor. Takes two two arguments. Creates
     *  component with given number and belonging to given domain.
     * @param n component number in particular domain. For instance, can represent
     * node number in particular domain.
     * @param d domain to which component belongs to */
    FEMComponent(int n, Domain *d) { number = n;
                                     domain = d; }
    /// virtual destructor
    virtual ~FEMComponent() { }                             // destructor

    /** Returns classType id of receiver. Intended for run time
     *  type checking. Every derived class have to overload this method.
     * @see cltypes.h include */
    virtual classType giveClassID() const { return FEMComponentClass; }
    /** Returns class name of the receiver.
     */
    virtual const char *giveClassName() const  = 0;
    /** Returns input record name of the receiver.
     */
    virtual const char *giveInputRecordName() const { return ( this->giveClassName() ); }
    /** Returns domain, which receiver belongs to.
     */
    Domain *giveDomain() const { return domain; }
    /// sets associated Domain
    virtual void         setDomain(Domain *d) { this->domain = d; }
    /** Returns component number of receiver.
     */
    int            giveNumber() const { return number; }
    /// sets number of receiver
    void setNumber(int _num) { this->number = _num; }
    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various fem components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflext updated numbering. The renumbering funciton is passed, which is supposed
     * to return an updated number of specified entyty type based on old number.
     */
    template< class T > void updateLocalNumbering( T *src, int ( T :: *renumberMethod )( int oldnum, EntityRenumberingScheme scheme ) ) { }
    /** Initializes receiver acording to object description stored in input record.
     *  This function is called immediately after creating object using
     * constructor. Input record can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /** Stores receiver state to output stream.
     *  Writes the FEMComponent class-id in order to allow test whether correct data are then restored.
     * @param stream output stream
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj special parameter, used only to send particular integration
     * point to material class version of this method. Except this
     * case, obj parameter is always NULL pointer.
     * @return contextIOResultType
     * @exception throws an ContextIOERR exception if error encountered
     */
    virtual contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /** Restores the receiver state previously written in stream.
     *  Readss the FEMComponent class-id in order to allow test consistency.
     * @see saveContext member function.
     * @return contextIOResultType
     * @exception throws an ContextIOERR exception if error encountered
     */
    virtual contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    // saves current context(state) into stream
    /** Allows programmer to test some internal data, before computation begins.
     *  For example, one may use this function, to ensure that element has material with
     * required capabilities is assigned to element. This must be done after all
     * mesh components are instanciated.
     * @return nonzero if receiver check is o.k. */
    virtual int    checkConsistency() { return 1; }
    /** Prints output of receiver to stream, for given time step */
    virtual void   printOutputAt(FILE *, TimeStep *) { }
    /** Interface requesting service */
    virtual Interface *giveInterface(InterfaceType) { return NULL; }
    /**@name error and warning reporting methods
     * These methods will print error (or warning) message using oofem default loggers.
     * Do not use these methods directly, to avoid specify file and line parameters.
     * More preferably, use these methods via corresponding OOFEM_CLASS_ERROR and OOFEM_CLASS_WARNING macros,
     * that will include file and line parameters automatically.
     *
     * Uses variable number of arguments, so a format string followed by optional argumens is expected
     * (according to printf conventions).
     * @param file  source file name, where error encountered (where error* function called)
     * @param line  source file line number, where error encountered
     */
    //@{
    /// prints error message and exits
    void error(const char *file, int line, const char *format, ...) const;
    /// prints warning message
    void warning(const char *file, int line, const char *format, ...) const;
    //@}
#ifdef __OOFEG
    virtual void   drawYourself(oofegGraphicContext &) { }
#endif
};

#endif // femcmpnn_h








