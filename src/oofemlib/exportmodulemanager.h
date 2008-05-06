/* $Header: /home/cvs/bp/oofem/oofemlib/src/exportmodulemanager.h,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   *********************************
//   *** CLASS ExportModuleManager ***
//   *********************************


#ifndef exportmodulemanager_h
#define exportmodulemanager_h

#include "alist.h"
#include "exportmodule.h"
#include "datareader.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <time.h>
#endif

class EngngModel;

/**
 * Class representing and implementing ExportModuleManager. It is attribute of EngngModel.
 * It manages the export output modules, which perform module - specific output oprations.
 */
class ExportModuleManager
{
private:
    /// module list
    AList< ExportModule > *moduleList;
    /// number of modules
    int numberOfModules;
    /// Associated Engineering model.
    EngngModel *emodel;

public:
    ExportModuleManager(EngngModel *emodel);
    ~ExportModuleManager();
    /**
     * Reads receiver description from input stream and creates corresponding modules components accordingly.
     * It scans input file, each line is assumed to be single record describing particular module.
     * The record line is converted to lowercase letters.
     * Corresponding component is created using ofType function.
     * After new output module object is created, its initializeForm member function is
     * called with its record as parameter.
     * @param inputStream input stream with domain description
     * @initString the e-model record containing export module manager record
     * @return nonzero if o.k.
     */
    int                instanciateYourself(DataReader *dr, InputRecord *ir);
    /**
     * Instanciates the receiver from input record. Called from instanciateYourself to initialize yourself
     * from corresponding record. Should be caled before instanciateYourself.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /**
     * Writes the output. Loops over all modules and calls corresponding doOutput module service.
     * @param tStep time step.
     */
    void              doOutput(TimeStep *tStep);
    /**
     * Initializes output manager. The corresponding initialize module services are called.
     */
    void              initialize();
    /**
     * Terminates the receiver, the corresponding terminate module services are called.
     */
    void              terminate();
    const char *giveClassName() const { return "IMLSolver"; }

protected:

    /**
     * Returns the required module.
     * @param num module number
     */
    ExportModule *giveExportModule(int num);
};


#endif // exportmodulemanager_h
