/* $Header: /home/cvs/bp/oofem/oofemlib/src/linearelasticmaterial.h,v 1.9 2003/04/06 14:08:25 bp Exp $ */
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


//   **************************************
//   *** CLASS LINEAR ELACSTIC MATERIAL ***
//   **************************************

#ifndef micromaterial_h
#define micromaterial_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "timestep.h"
#include "oofem_limits.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "datastream.h"
#include "contextioerr.h"

class MicroMaterialStatus : public StructuralMaterialStatus {

public:
/// constructor
    MicroMaterialStatus(int, Domain *d, GaussPoint *gp);

/// destructor
~MicroMaterialStatus();
void  initTempStatus();
void  updateYourself(TimeStep *atTime);
void  printOutputAt(FILE *file, TimeStep *tStep);

const char *giveClassName() const { return "MicroMaterialStatus"; }
classType  giveClassID() const  { return MicroMaterialStatusClass; }


contextIOResultType  saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    
contextIOResultType  restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

protected:
};


/**
 * This class is an base class for microproblem. The microproblem represents itself a problem which is solved separately from the macroproblem with appropriate boundary conditions. Stiffness matrix of microproblem is condensed to provide stiffness matrix for macroelement.
Also general implementation of giveRealStressVector service is provided,
computing the stress increment vector from strain increment multiplied by
stiffness.
 */
class MicroMaterial : public StructuralMaterial
{
public:
    /// Constructor
    MicroMaterial(int n, Domain *d);
    /// Destructor
    ~MicroMaterial() { };

    char inputFileNameMicro [ MAX_FILENAME_LENGTH ];

    IRResultType initializeFrom (InputRecord* ir);

    //int MicroMaterial :: hasMaterialModeCapability(MaterialMode mode);

const char *giveClassName() const { return "MicroMaterial"; }


void giveRealStressVector (FloatArray& answer,  MatResponseForm, GaussPoint*, const FloatArray&, TimeStep*);


MaterialStatus *CreateStatus(GaussPoint *gp) const;

///pointer to the underlying micro problem
EngngModel *problemMicro;

protected:


};


#endif // micromaterial_h
