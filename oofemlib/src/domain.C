/* $Header: /home/cvs/bp/oofem/oofemlib/src/domain.C,v 1.31.4.2 2004/05/14 13:45:27 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

//   file DOMAIN.CC

#include "domain.h"
#include "element.h"
#include "timestep.h"
#include "node.h"
#include "elementside.h"
#include "material.h"
#include "crosssection.h"
//#include "yieldcriteria.h"
#include "load.h"
#include "initial.h"
#include "loadtime.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
//#include <clock.h>
#endif
#include "clock.h"
#include "debug.h"
#include "verbose.h"
#include "strreader.h"
#include "cltypes.h"
#include "conTable.h"
#include "outputmanager.h"
#include "dummylocalizer.h"
#include "octreelocalizer.h"
#include "datareader.h"
#include "util.h"
#include "nodalrecoverymodel.h"
#include "nonlocalbarrier.h"
#include "usrdefsub.h"
#include "logger.h"

#ifdef __PARALLEL_MODE
#include "parallel.h"
#include "processcomm.h"
#include "loadballancer.h"
#include "datastream.h"
#include "communicator.h"
#include "parmetisloadballancer.h"
#endif

#ifndef __MAKEDEPEND
#include <string.h>
#include <stdarg.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#include <ctype.h>
#endif

Domain :: Domain (int n, int serNum, EngngModel* pm) : defaultNodeDofIDArry(), defaultSideDofIDArry()
   // Constructor. Creates a new domain.
{
 this->engineeringModel = pm;
 this->number   = n;
 this->serialNumber = serNum;

 elementList           = new AList<Element>(0) ;
 dofManagerList        = new AList<DofManager>(0) ;
 materialList          = new AList<Material>(0) ;
 bcList                = new AList<GeneralBoundaryCondition>(0) ;
 icList                = new AList<InitialCondition>(0) ;
 loadTimeFunctionList  = new AList<LoadTimeFunction>(0) ;
 crossSectionList      = new AList<CrossSection>(0) ;
 nonlocalBarierList    = new AList<NonlocalBarrier>(0) ;
 // yieldCriteriaList     = new AList(0) ;

 numberOfDefaultDofsPerNode = -1 ;
 numberOfDefaultDofsPerSide = -1 ;
 dType                 = _unknownMode ;
 
 connectivityTable     = NULL;
 spatialLocalizer      = NULL;
 outputManager         = new OutputManager (this);
 smoother = NULL;

 nonlocalUpdateStateCounter = 0;
#ifdef __PARALLEL_MODE
 dmanMapInitialized = false;
#endif
}

Domain :: ~Domain ()
   // Destructor.
{
   delete elementList ;
   delete dofManagerList ;
   delete materialList ;
   delete bcList ;
   delete icList ;
   delete loadTimeFunctionList ;
   delete crossSectionList ;
   delete nonlocalBarierList ;

   delete connectivityTable;
   delete spatialLocalizer;
   delete outputManager;
   if (smoother) delete smoother;

}


Element*  Domain :: giveElement (int n)
   // Returns the n-th element. Generates error if it is not defined yet.
{

   if (elementList -> includes(n))
      return elementList->at(n) ;
   else {
     _error2 ("giveElement: undefined element (%d)",n);
     // elem = (Element*) Element(n,this).typed() ;
     // elementList -> put(n,elem) ;}
   }
   return NULL ;
}

/*
FILE*  Domain :: giveInputStream ()
   // Returns an input stream on the data file of the receiver.
{
   if (inputStream)
      return inputStream ;
   else {
   fprintf (stderr,"\nDomain->giveInputStream: Internal error\a\n") ;
   exit (1);}
  
  return inputStream ;
}

FILE*  Domain :: giveOutputStream ()
   // Returns an output stream on the data file of the receiver.
{
   if (! outputStream) {
   fprintf (stderr,"\nDomain->giveOutputStream: Internal error\a\n") ;
   exit (1);
  }
   return outputStream ;
}
*/
Load*  Domain :: giveLoad (int n)
   // Returns the n-th load. Generates the error if not defined.
{
  Load* answer;

  if (bcList -> includes(n)) {
      answer = dynamic_cast<Load*> (bcList->at(n)) ;
      if (answer) {
        return answer;
      } else {
        _error2 ("giveLoad: cannot cast boundary condition %d to Load class", n);
      }
   } else {
     _error2 ("giveLoad: undefined load (%d)",n);
//      load = (Load*) Load(n,this).typed() ;
//      loadList -> put(n,load) ;}
   }
   return NULL ;
}

GeneralBoundaryCondition*  Domain :: giveBc (int n)
   // Returns the n-th bc. Generates the error if not defined.
{
  if (bcList -> includes(n)) {
    return bcList->at(n) ;
  } else {
    _error2 ("giveBc: undefined bc (%d)",n);
  }
  return NULL ;
}

InitialCondition*  Domain :: giveIc (int n)
   // Returns the n-th ic. Generates the error if not defined.
{
  if (icList -> includes(n)) {
    return icList->at(n) ;
  } else {
    _error2 ("giveIc: undefined ic (%d)",n);
  }
  return NULL ;
}


LoadTimeFunction*  Domain :: giveLoadTimeFunction (int n)
   // Returns the n-th load-time function. Creates this fuction if it does
   // not exist yet.
{
   if (loadTimeFunctionList -> includes(n))
      return loadTimeFunctionList->at(n) ;
   else {
     _error2 ("giveLoadTimeFunction: undefined load-time function (%d)",n);
     //      ltf = (LoadTimeFunction*) LoadTimeFunction(n,this).typed() ;
     //      loadTimeFunctionList -> put(n,ltf) ;}
   }
   return NULL ;
}


Material*  Domain :: giveMaterial (int n)
   // Returns the n-th material. Creates this material if it does not exist
   // yet.
{

   if (materialList -> includes(n))
      return materialList -> at(n) ;
   else {
     _error2 ("giveMaterial: undefined material (%d)",n);
     //      mat = new Material(n,this) ;
     //      materialList  -> put(n,mat) ;}
   }
   return NULL ;
}


Node*  Domain :: giveNode (int n)
   // Returns the n-th node. Creates this node if it does not exist yet.
{
   DofManager *node=NULL ;

   if (dofManagerList -> includes(n)) {
      node = dofManagerList -> at(n) ;
   if ((node->giveClassID() != NodeClass) && (node->giveClassID() != RigidArmNodeClass) && (node->giveClassID() != HangingNodeClass)) {
    _error2 ("giveNode: incompatible type of dofManager %d, can not convert",n);
   }
  } else {
   _error2 ("giveNode: undefined dofManager (%d)",n);
   //      node = new Node(n,this) ;
   //      nodeList  -> put(n,node) ;}
  }
   return (Node*) node ;
 }

ElementSide*  Domain :: giveSide (int n)
   // Returns the n-th element side.
{
   DofManager *side = NULL ;

   if (dofManagerList -> includes(n)) {
      side = dofManagerList -> at(n) ;
   if (side->giveClassID() != ElementSideClass) {
    _error2 ("giveSide: incompatible type of dofManager %d, can not convert",n);
   }
  } else {
   _error2 ("giveSide: undefined dofManager (%d)",n);
  }
   return (ElementSide*) side ;
 }

DofManager*  Domain :: giveDofManager (int n)
   // Returns the n-th node. Creates this node if it does not exist yet.
{
   if (dofManagerList -> includes(n)) {
   return dofManagerList -> at(n) ;
  } else {
   _error2 ("giveDofManager: undefined dofManager (%d)",n);
   //      node = new Node(n,this) ;
   //      nodeList  -> put(n,node) ;}
  }
   return NULL ;
 }



CrossSection*  Domain :: giveCrossSection (int n)
   // Returns the n-th cross section.
   // yet.
{

   if (crossSectionList -> includes(n))
      return crossSectionList -> at(n) ;
   else {
     _error2 ("giveCrossSection: undefined cross section (%d)",n);
   }
   return NULL;
}


NonlocalBarrier*  Domain :: giveNonlocalBarrier (int n)
   // Returns the n-th NonlocalBarrier.
{

   if (nonlocalBarierList -> includes(n))
      return nonlocalBarierList -> at(n) ;
   else {
     _error2 ("giveNonlocalBarrier: undefined barrier (%d)",n);
   }
   return NULL;
}



/*
YieldCriteria*  Domain :: giveYieldCriteria (int n)
   // Returns the n-th yieldCriteria
   // yet.
{
   YieldCriteria* yieldCriteria ;

   if (yieldCriteriaList -> includes(n))
      yieldCriteria = (YieldCriteria*) yieldCriteriaList -> at(n) ;
   else {
     _errori ("giveYieldCriteria: No such cross section defined: ",n);
   }
   return yieldCriteria ;
}
*/

EngngModel*  Domain :: giveEngngModel ()
   // Returns the time integration algorithm. Creates it if it does not
   // exist yet.
{
 if (engineeringModel)
      return  engineeringModel;
   else {
  _error("giveEngngModel: Not defined");
   }
 
 return NULL;  
}


int  Domain :: instanciateYourself (DataReader* dr)
   // Creates all objects mentioned in the data file.

{
 const char *__keyword, *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                               // Required by IR_GIVE_FIELD macro

   int i, num ;
   char name [MAX_NAME_LENGTH];
   int nnode,nelem,nmat,nload,nic,nloadtimefunc,ncrossSections,nbarrier ;
   DofManager *node ;
   Element *elem ;
   Material *mat ;
   GeneralBoundaryCondition *load;
   InitialCondition *ic;
   LoadTimeFunction *ltf;
  CrossSection *crossSection;
  NonlocalBarrier *barrier;
  FILE *outputStream = this->giveEngngModel()->giveOutputStream();

 // read type of Domain to be solved
 InputRecord* ir = dr->giveInputRecord (DataReader::IR_domainRec, 1);
 __keyword = "domain"; result = ir->giveField(name, MAX_NAME_LENGTH, IFT_Domain_type, __keyword);
 if (result != IRRT_OK) IR_IOERR (giveClassName(), __proc, IFT_Domain_type, __keyword, ir, result);
 ir->finish();
 
#  ifdef VERBOSE
 VERBOSE_PRINT0("Instanciating domain ", this->number);
#  endif
 
 resolveDomainDofsDefaults (name);
 fprintf (outputStream,"Domain type: %s, default ndofs per node is %d, per side is %d\n\n\n",
      name, giveNumberOfDefaultNodeDofs (),giveNumberOfDefaultSideDofs() );
 
 // read otput manager record
 ir = dr->giveInputRecord (DataReader::IR_outManRec, 1);
 outputManager->initializeFrom (ir);
 ir->finish();

 // read domain description
 ir = dr->giveInputRecord (DataReader::IR_domainCompRec, 1);
 IR_GIVE_FIELD (ir, nnode, IFT_Domain_ndofman, "ndofman"); // Macro
 IR_GIVE_FIELD (ir, nelem, IFT_Domain_nelem, "nelem"); // Macro
 IR_GIVE_FIELD (ir, ncrossSections, IFT_Domain_ncrosssect, "ncrosssect"); // Macro
 IR_GIVE_FIELD (ir, nmat, IFT_Domain_nmat, "nmat"); // Macro
 IR_GIVE_FIELD (ir, nload, IFT_Domain_nbc, "nbc"); // Macro
 IR_GIVE_FIELD (ir, nic, IFT_Domain_nic, "nic"); // Macro
 IR_GIVE_FIELD (ir, nloadtimefunc, IFT_Domain_nloadtimefunct, "nltf"); // Macro

 // read optional number of nonlocalBarriers
 nbarrier = 0; 
 __keyword = "nbarrier"; result = ir->giveOptionalField(nbarrier,  IFT_Domain_nbarrier, __keyword);

 // read nodes
 dofManagerList -> growTo(nnode) ;
 for (i=0; i < nnode ; i++) {
   ir = dr->giveInputRecord (DataReader::IR_dofmanRec, i+1);
   // read type of dofManager
   //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
   IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

   (node = (DofManager *) 
    (DofManager(num,this).ofType(name)))->initializeFrom(ir);

   // check number
   if ((num < 1) || (num > nnode)) 
     _error2 ("instanciateYourself: Invalid dofManager number (num=%d)", num);
   if (!dofManagerList->includes(num)) {
     dofManagerList->put(num,node) ;
   } else {
     _error2 ("instanciateYourself: Dofmanager entry already exist (num=%d)", num);
   }

   //dofManagerList->put(i+1,node) ;
   ir->finish();
 }
  
#  ifdef VERBOSE
  VERBOSE_PRINT0("Instanciated nodes & sides ",nnode)
#  endif

// read elements
   elementList -> growTo(nelem) ;
   for (i=0; i < nelem ; i++) { 
     ir = dr->giveInputRecord (DataReader::IR_elemRec, i+1);
     // read type of element
     //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
     IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

     (elem = (Element *) 
      (Element(num,this).ofType(name)))->initializeFrom(ir);

     // check number
     if ((num < 1) || (num > nelem)) 
       _error2 ("instanciateYourself: Invalid element number (num=%d)", num);
     if (!elementList->includes(num)) {
       elementList->put(num,elem) ;
     } else {
       _error2 ("instanciateYourself: element entry already exist (num=%d)", num);
     }

   //elementList->put(i+1,elem) ;
   ir->finish();
   }

#  ifdef VERBOSE
 VERBOSE_PRINT0("Instanciated elements ",nelem);
#  endif

// read cross sections
 crossSectionList -> growTo(ncrossSections);
   for (i=0; i < ncrossSections ; i++) {
     ir = dr->giveInputRecord (DataReader::IR_crosssectRec, i+1);
     //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
     IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

     (crossSection  = (CrossSection *) 
      (CrossSection(num,this).ofType(name)))->initializeFrom(ir);

     // check number
     if ((num < 1) || (num > ncrossSections)) 
       _error2 ("instanciateYourself: Invalid crossSection number (num=%d)", num);
     if (!crossSectionList->includes(num)) {
       crossSectionList->put(num,crossSection) ;
     } else {
       _error2 ("instanciateYourself: crossSection entry already exist (num=%d)", num);
     }
     //crossSectionList->put(i+1,crossSection) ;
     ir->finish();
   }   

#  ifdef VERBOSE
 VERBOSE_PRINT0 ("Instanciated cross sections ",ncrossSections) 
#  endif

  // read materials
  materialList -> growTo (nmat);
  for (i=0; i < nmat ; i++) { 
    ir = dr->giveInputRecord (DataReader::IR_matRec, i+1);
    // read type of material
    //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
    IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

    (mat  = (Material *) 
     (Material(num,this).ofType(name)))->initializeFrom(ir);

    // check number
    if ((num < 1) || (num > nmat)) 
      _error2 ("instanciateYourself: Invalid material number (num=%d)", num);
    if (!materialList->includes(num)) {
      materialList->put(num, mat) ;
    } else {
      _error2 ("instanciateYourself: material entry already exist (num=%d)", num);
    }
    //materialList->put(i+1,mat) ;
    ir->finish();
  }   

#  ifdef VERBOSE
 VERBOSE_PRINT0 ("Instanciated materials ",nmat) 
#  endif

  if (nbarrier) {
    // read barriers
    nonlocalBarierList -> growTo (nbarrier);
    for (i=0; i < nbarrier ; i++)
      { 
        ir = dr->giveInputRecord (DataReader::IR_nlocBarRec, i+1);
        // read type of load
        //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

        barrier = CreateUsrDefNonlocalBarrierOfType (name, num, this);
        barrier->initializeFrom(ir);

        // check number
        if ((num < 1) || (num > nbarrier)) 
          _error2 ("instanciateYourself: Invalid barrier number (num=%d)", num);
        if (!nonlocalBarierList->includes(num)) {
          nonlocalBarierList->put(num, barrier) ;
        } else {
          _error2 ("instanciateYourself: barrier entry already exist (num=%d)", num);
        }
        //nonlocalBarierList->put(i+1,barrier) ;
        ir->finish();
      }
#  ifdef VERBOSE
    VERBOSE_PRINT0 ("Instanciated barriers ",nbarrier) ;
#  endif
   
  }


// read boundary conditions
 bcList -> growTo (nload);
 for (i=0; i < nload ; i++) {  
   ir = dr->giveInputRecord (DataReader::IR_bcRec, i+1);
   // read type of load
   //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
   IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

   (load = ( GeneralBoundaryCondition*) 
    (GeneralBoundaryCondition(num,this).ofType(name)))->initializeFrom(ir);

   // check number
   if ((num < 1) || (num > nload)) 
     _error2 ("instanciateYourself: Invalid boundary condition number (num=%d)", num);
   if (!bcList->includes(num)) {
     bcList->put(num, load) ;
   } else {
     _error2 ("instanciateYourself: boundary condition entry already exist (num=%d)", num);
   }
   //loadList->put(i+1,load) ;
   ir->finish();
 }

#  ifdef VERBOSE
 VERBOSE_PRINT0 ("Instanciated BCs ",nload)
#  endif

// read initial conditions
 icList -> growTo (nic);
 for (i=0; i < nic ; i++) {  
   ir = dr->giveInputRecord (DataReader::IR_icRec, i+1);
   // read type of load
   //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
   IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

   ic = new InitialCondition(num,this);
   if (ic) {
     ic->initializeFrom(ir);
   } else _error2 ("instanciateYourself: Creation of IC no. %d failed", num);

   // check number
   if ((num < 1) || (num > nic)) 
     _error2 ("instanciateYourself: Invalid initial condition number (num=%d)", num);
   if (!icList->includes(num)) {
     icList->put(num, ic) ;
   } else {
     _error2 ("instanciateYourself: initial condition entry already exist (num=%d)", num);
   }
   ir->finish();
 }

#  ifdef VERBOSE
 VERBOSE_PRINT0 ("Instanciated ICs ",nic)
#  endif


 // read load time functions
 loadTimeFunctionList -> growTo (nloadtimefunc);
 for (i=0; i < nloadtimefunc ; i++) {  
   ir = dr->giveInputRecord (DataReader::IR_ltfRec, i+1);
   // read type of ltf
   //__keyword = NULL; result = ir->giveField(name, MAX_NAME_LENGTH, __keyword);
   IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num, MAX_NAME_LENGTH);

   (ltf  = (LoadTimeFunction *) 
    (LoadTimeFunction(num,this).ofType(name)))->initializeFrom(ir);

   // check number
   if ((num < 1) || (num > nloadtimefunc)) 
     _error2 ("instanciateYourself: Invalid LoadTimeFunction number (num=%d)", num);
   if (!loadTimeFunctionList->includes(num)) {
     loadTimeFunctionList->put(num, ltf) ;
   } else {
     _error2 ("instanciateYourself: LoadTimeFunction entry already exist (num=%d)", num);
   }
   //loadTimeFunctionList->put(i+1,ltf) ;
   ir->finish();
 }   

#  ifdef VERBOSE
 VERBOSE_PRINT0 ("Instanciated load-time fncts ",nloadtimefunc)
#  endif

  return 1;
}


void  Domain :: error (char* file, int line, char *format, ...) 
{
  char buffer[MAX_ERROR_MSG_LENGTH];
	va_list args;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
  
  __OOFEM_ERROR3 (file, line, "Class: Domain, number: %d\n%s",number,buffer);
}


void  Domain :: warning (char* file, int line, char *format, ...) 
{
  char buffer[MAX_ERROR_MSG_LENGTH];
	va_list args;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
  
  __OOFEM_WARNING3 (file, line, "Class: Domain, number: %d\n%s",number,buffer);
}


const IntArray&
Domain :: giveDefaultNodeDofIDArry ()
{
// returns default DofID array, defining physical meaning of partucular DOFs
// in Node Dof collection
 if (this->defaultNodeDofIDArry.giveSize()) {
  return defaultNodeDofIDArry;
 }


  if(dType == _2dPlaneStressRotMode) {
  defaultNodeDofIDArry.resize(3);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_v;defaultNodeDofIDArry.at(3)=R_w;
 }
  else if(dType == _2dPlaneStressMode) {
  defaultNodeDofIDArry.resize (2);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_v;
 }
  else if(dType == _PlaneStrainMode) {
  defaultNodeDofIDArry.resize (2);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_v;
 }
  else if  (dType == _3dMode) {
  defaultNodeDofIDArry.resize (3);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_v; defaultNodeDofIDArry.at(3)=D_w;
 }
 else if (dType == _3dAxisymmMode) {
  defaultNodeDofIDArry.resize (3);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_v; defaultNodeDofIDArry.at(3)=R_w;
 }
  else if  (dType == _2dMindlinPlateMode) {
  defaultNodeDofIDArry.resize (3);
    defaultNodeDofIDArry.at(1)=D_w; defaultNodeDofIDArry.at(2)=R_u; defaultNodeDofIDArry.at(3)=R_v;
 }
 else if ( dType == _3dShellMode) {
  defaultNodeDofIDArry.resize (6);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_v; defaultNodeDofIDArry.at(3)=D_w;
  defaultNodeDofIDArry.at(4)=R_u; defaultNodeDofIDArry.at(5)=R_v; defaultNodeDofIDArry.at(6)=R_w;
 }
  else if  (dType == _2dTrussMode) {
  defaultNodeDofIDArry.resize (2);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_w;
 }
  else if  (dType == _1dTrussMode) {
  defaultNodeDofIDArry.resize(1);
    defaultNodeDofIDArry.at(1)=D_u; 
 }
  else if  (dType == _2dBeamMode) {
  defaultNodeDofIDArry.resize (3);
    defaultNodeDofIDArry.at(1)=D_u; defaultNodeDofIDArry.at(2)=D_w; defaultNodeDofIDArry.at(3)=R_v;
 }
  else if  (dType == _HeatTransferMode) {
  defaultNodeDofIDArry.resize (1);
    defaultNodeDofIDArry.at(1)=T_f;
 }  
  else if  (dType == _HeatMass1Mode) {
  defaultNodeDofIDArry.resize (2);
    defaultNodeDofIDArry.at(1)=T_f;
    defaultNodeDofIDArry.at(2)=C_1;
  }  else if (dType == _2dIncompressibleFlow) {
    defaultNodeDofIDArry.resize (3);
    defaultNodeDofIDArry.at(1)=V_u;
    defaultNodeDofIDArry.at(2)=V_v;
    defaultNodeDofIDArry.at(3)=P_f;
  } else {
    _error("giveDefaultNodeDofIDArry : unknown domainType");
  }
 return defaultNodeDofIDArry;
}



int Domain ::  giveNumberOfDefaultNodeDofs () 
//
// returns default number of dofs for one node
// this number depend on type of problem (2dplane-stress, 3d truss, 3d, 2d beam etc.)
// returns member data  numberOfDefaultDofsPerNode.
// numberOfDefaultDofsPerNode is initialized in initiazeFrom subroutine.
// 
{
  if ( numberOfDefaultDofsPerNode == -1) {
    OOFEM_LOG_WARNING("Domain ::  giveNumberOfDefaultNodeDofs : Number of Default Dofs per Node is not specified, using default 6 instead\a\n");
    return (numberOfDefaultDofsPerNode = 6);
  } else return  numberOfDefaultDofsPerNode ;
}


const IntArray&
Domain :: giveDefaultSideDofIDArry ()
{
// returns default DofID array, defining physical meaning of partucular DOFs
// in side Dof collection

  // IntArray* answer;

/*
  if(dType == _2dPlaneStressRotMode) {
  answer = new IntArray (3);
    answer->at(1)=D_u; answer->at(2)=D_v;answer->at(3)=R_w;
 }
  else if(dType == _2dPlaneStressMode) {
  answer = new IntArray (2);
    answer->at(1)=D_u; answer->at(2)=D_v;
 }
  else if  (dType == _3dMode) {
  answer = new IntArray (3);
    answer->at(1)=D_u; answer->at(2)=D_v; answer->at(3)=D_w;
 }
 else if (dType == _3dAxisymmMode) {
  answer = new IntArray (3);
    answer->at(1)=D_u; answer->at(2)=D_v; answer->at(3)=R_w;
 }
  else if  (dType == _2dMindlinPlateMode) {
  answer = new IntArray (3);
    answer->at(1)=D_w; answer->at(2)=R_u; answer->at(3)=R_v;
 }
 else if ( dType == _3dShellMode) {
  answer = new IntArray (5);
    answer->at(1)=D_u; answer->at(2)=D_v; answer->at(3)=D_w;
  answer->at(4)=R_u; answer->at(5)=R_v;
 }
  else if  (dType == _2dTrussMode) {
  answer = new IntArray (2);
    answer->at(1)=D_u; answer->at(2)=D_w;
 }
  else if  (dType == _1dTrussMode) {
  answer = new IntArray (1);
    answer->at(1)=D_u; 
 }
  else if  (dType == _2dBeamMode) {
  answer = new IntArray (3);
    answer->at(1)=D_u; answer->at(2)=D_w; answer->at(3)=R_v;
 }
  else if  (dType == _2dHeatMode) {
  answer = new IntArray (1);
    answer->at(1)=T_f;
 }  
  else {
    _error("Domain : Domain type name of unknown type\a\n");
    return NULL;
  }
 return answer;
 */

  _error("giveDefaultSideDofIDArry : unknown domainType");
  defaultSideDofIDArry.resize(0);
  return defaultSideDofIDArry;
}



int Domain ::  giveNumberOfDefaultSideDofs () 
//
// returns default number of dofs for one side
// this number depend on type of problem (2dplane-stress, 3d truss, 3d, 2d beam etc.)
// returns member data  numberOfDefaultDofsPerNode.
// numberOfDefaultDofsPerNode is initialized in initializeFrom subroutine.
// 
{
  if ( numberOfDefaultDofsPerSide == -1) {
    _warning("giveNumberOfDefaultSideDofs: Number of Default Dofs per Side is not specified, using default 0 instead");
    return (numberOfDefaultDofsPerSide = 0);
  } else return  numberOfDefaultDofsPerSide ;
}

 



void Domain ::  resolveDomainDofsDefaults(char* typeName)
//
// resolves default number of dofs per node according to domain type name.
// and also resolves default dof mask according to domain type.
//
{
 numberOfDefaultDofsPerSide = 0;

  if(!strncasecmp(typeName,"2dplanestressrot",16)) {
    dType = _2dPlaneStressRotMode; 
  numberOfDefaultDofsPerNode = 3;
 }
  else if(!strncasecmp(typeName,"2dplanestress",12)) {
    dType = _2dPlaneStressMode; 
  numberOfDefaultDofsPerNode = 2;
 }
  else if(!strncasecmp(typeName,"planestrain",11)) {
    dType = _PlaneStrainMode; 
  numberOfDefaultDofsPerNode = 2;
 }
 else if (! strncasecmp(typeName,"3daxisymm",9)) {
    dType = _3dAxisymmMode;
  numberOfDefaultDofsPerNode = 3;
 }
  else if  (! strncasecmp(typeName,"2dmindlinplate",14)) {
    dType = _2dMindlinPlateMode;
  numberOfDefaultDofsPerNode = 3;
 }
 else if (! strncasecmp(typeName,"3dshell",7)) {
    dType = _3dShellMode;
  numberOfDefaultDofsPerNode = 6;
 }
  else if  (! strncasecmp(typeName,"2dtruss",7)){
    dType = _2dTrussMode; 
  numberOfDefaultDofsPerNode = 2;
 }
  else if  (! strncasecmp(typeName,"1dtruss",7)){
    dType = _1dTrussMode; 
  numberOfDefaultDofsPerNode = 1;
 }
  else if  (! strncasecmp(typeName,"2dbeam",6)) {
    dType = _2dBeamMode;  
  numberOfDefaultDofsPerNode = 3;
 }
  else if  (! strncasecmp(typeName,"3d",2))     {
    dType = _3dMode;
  numberOfDefaultDofsPerNode = 3;
 }
  else if  (! strncasecmp(typeName,"heattransfer",11)) {
    dType = _HeatTransferMode;
  numberOfDefaultDofsPerNode = 1;
 }  
  else if  (! strncasecmp(typeName,"hema1",5)) {
    dType = _HeatMass1Mode;
  numberOfDefaultDofsPerNode = 2;
  } else if (! strncasecmp(typeName,"2dincompflow",12)) {
    dType = _2dIncompressibleFlow;
    numberOfDefaultDofsPerNode = 3;
  } else { 
    _error("resolveDomainDofsDefaults : unknown domainType");
    return;
  }
}
    

#ifdef __OOFEG

void Domain :: drawYourself (oofegGraphicContext& context) 
//
// shows graphics representation of domain, respecting mode
//
{
 
 OGC_PlotModeType plotMode = context.giveIntVarPlotMode();
 if ((plotMode == OGC_nodeAnnotation) || (plotMode == OGC_nodeGeometry) || (plotMode == OGC_essentialBC) ||
     (plotMode == OGC_naturalBC) || (plotMode == OGC_nodeScalarPlot) || (plotMode == OGC_nodeVectorPlot))
    this->drawNodes (context);
 else 
  this->drawElements (context);
}

void Domain :: drawElements (oofegGraphicContext& context) {
//
// steps through element array and calls element(i)->show(mode,this);
//
  for (int i=1; i <= this->giveNumberOfElements() ; i++)
    {  
      this->giveElement(i)->drawYourself(context);
    }
}

void  Domain :: drawNodes (oofegGraphicContext& context) {
//
// steps through element array and calls element(i)->show(mode,this);
//
  int nnodes = this->giveNumberOfDofManagers();
  for (int i=1; i <= nnodes ; i++)
    {  
      this->giveDofManager(i)->drawYourself(context);
    }
}

#endif


NodalRecoveryModel* 
Domain :: giveSmoother ()
{
 return this->smoother;
}

void
Domain :: setSmoother (NodalRecoveryModel* smoother, int destroyOld)
{
 if (destroyOld && this->smoother) delete this->smoother;
 this->smoother = smoother;
}




ConnectivityTable * Domain :: giveConnectivityTable () 
// 
// return connectivity Table - if no defined - creates new one
//
{
  if (connectivityTable == NULL) connectivityTable = new ConnectivityTable(this);
  return connectivityTable;
}


SpatialLocalizer * Domain :: giveSpatialLocalizer () 
// 
// return connectivity Table - if no defined - creates new one
//
{

//  if (spatialLocalizer == NULL) spatialLocalizer = new DummySpatialLocalizer(1, this);
 if (spatialLocalizer == NULL) spatialLocalizer = new OctreeSpatialLocalizer(1, this);
  return spatialLocalizer;
}




int Domain ::  giveCorrespondingCoordinateIndex (int idof)
//
// find corresponding coordinate axis to idof
// if no - coordinate axis corespond to idof returns 0;
// 
// if idof corresponds to displacement in direction of axis i then finction returns i
// otherwise 0;
//
{
  switch (dType) {
  case _2dBeamMode:
    if (idof == 1) return 1; else if(idof == 2) return 3;
    return 0;

  case _2dPlaneStressMode:
 case _PlaneStrainMode:
    if (idof == 1) return 1; else if(idof == 2) return 2;
    return 0;

  case _2dPlaneStressRotMode:
    if (idof == 1) return 1; else if(idof == 2) return 2;
    return 0;

  case _2dTrussMode:
    if (idof == 1) return 1; else if(idof == 2) return 3;
    return 0;

  case _1dTrussMode:
    if (idof == 1) return 1; 
    return 0;

  case _2dMindlinPlateMode:
    if (idof == 1) return 3; 
    return 0;

 case _3dMode:
  if (idof == 1) return 1; else if (idof == 2) return 2;
  else if (idof == 3) return 3;
  return 0;

 case _3dAxisymmMode:
  if (idof == 1) return 1; else if (idof == 2) return 2;
  return 0;

 case _3dShellMode:
  if (idof == 1) return 1; else if (idof == 2) return 2;
  else if (idof == 3) return 3;
  return 0;


  
  default:
    _error ("giveCorrespondingCoordinateIndex : unsupported domain type");
  }
return 0;
}

/*
Domain :: giveCorrespondingDofID (int idof)
{
// returns corresponding DofId to idof-th dof in node
// respecting current domain mode.
// if no corresponding dofID exists returns (Err_dof = 0)
//
 
 switch (dType) {
 case _2dBeamMode:
  if     (idof == 1) return D_u; 
  else if(idof == 2) return D_w;
  else if(idof == 3) return R_v;
  break ;
 case _2dPlaneStressMode:
  if     (idof == 1) return D_u;
  else if(idof == 2) return D_v;
  break;
 case _2dTrussMode:
  if     (idof == 1) return D_u;
  else if(idof == 2) return D_v;
  break;
 case _1dTrussMode:
  if     (idof == 1) return D_u;
  break;
 case _2dMindlinPlateMode:
  if     (idof == 1) return D_w;
  else if(idof == 2) return R_u;
  else if(idof == 3) return R_v;
  break;
 case _3dMode:
  if     (idof == 1) return D_u;
  else if(idof == 2) return D_v;
  else if(idof == 3) return D_w;
  break;
 case _2dHeatMode:
   if     (idof == 1) return T_f;
   break;
 default:
  _error ("giveCorrespondingDofID : udefined iDof for selected domainType");
 }
 return Err_dof;
 
}
*/

int
Domain :: checkConsistency ()
//
// checks internal consistency
// 
// many parameters are checked at run-time during computation
//
// this function transverse tree of all objects and invokes
// checkConsistency on this objects
// currently this function checks noly consistency
// of internal object structures, mainly whether referenced other objects
// are having required support
// 
{
  int i, result=1;
  int nnode, nelem, nmat;
  
  nnode = this-> giveNumberOfDofManagers ();
  nelem = this-> giveNumberOfElements();
  nmat  = this-> giveNumberOfMaterialModels();

  for (i=1; i <= nnode ; i++)
  result &= this->giveDofManager(i)->checkConsistency();

  for (i=1; i <= nelem ; i++)
  result &= this->giveElement(i)->checkConsistency();

  for (i=1; i <= nmat ; i++)
  result &= this->giveMaterial(i)->checkConsistency();

  return result;
}

ErrorEstimator*
Domain::giveErrorEstimator () {
 return engineeringModel->giveDomainErrorEstimator (this->number);
}


#ifdef __PARALLEL_MODE

/***********************************
 Load Ballancing methods
***********************************/

/* packs data for remote partition
   please note, that if entity has to be physically moved (not updated) then
   all links (such as element node numbers) should be packed in global numbering !!
   the remote local numbering could not be recovered on local partition!
 */
int
Domain::packMigratingData (LoadBallancer* lb, ProcessCommunicator& pc) 
{
  int myrank = this->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int idofman, ndofman;
  classType dtype;
  DofManager* dofman;
  LoadBallancer::DofManMode dmode;

  //  **************************************************
  //  Pack migrating data to remote partition 
  //  **************************************************

  // pack dofManagers
  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);
  // loop over dofManagers
  ndofman = this->giveNumberOfDofManagers();
  for (idofman=1; idofman <= ndofman; idofman++) {
    dofman = this->giveDofManager (idofman);
    dmode = lb->giveDofManState(idofman);
    dtype = dofman->giveClassID();
    // sync data to remote partition 
    // if dofman already present on remote partition then there is no need to sync
    //if ((lb->giveDofManPartitions(idofman)->findFirstIndexOf(iproc))) {
    if ((lb->giveDofManPartitions(idofman)->findFirstIndexOf(iproc)) &&
       (!dofman->givePartitionList()->findFirstIndexOf(iproc))) { 
      pcbuff->packInt (dtype);
      pcbuff->packInt (dmode);
      pcbuff->packInt (dofman->giveGlobalNumber());

      // pack dofman state (this is the local dofman, not available on remote)
      /* this is a potential performance leak, sending shared dofman to a partition,
         in which is already shared does not require to send context (is already there)
         here for simplicity it is always send */
      dofman->saveContext (&pcDataStream, CM_Definition | CM_State | CM_UnknownDictState);
      // send list of new partitions
      pcbuff->packIntArray (*(lb->giveDofManPartitions(idofman)));
    }
  }

  // pack end-of-dofman-section record
  pcbuff->packInt (LOADBALLANCER_END_DATA);
  
  int ielem, nelem = this->giveNumberOfElements();
  
  Element* elem;

  for (ielem=1; ielem<=nelem; ielem++) { // begin loop over elements
    elem = this->giveElement(ielem);
    if ((elem->giveParallelMode() == Element_local) && 
	(lb->giveElementPartition(ielem) == iproc)) {
      // pack local element (node numbers shuld be global ones!!!)
      // pack type
      pcbuff->packInt (elem->giveClassID());
      // nodal numbers shuld be packed as global !!
      elem->saveContext (&pcDataStream, CM_Definition | CM_DefinitionGlobal | CM_State);
    }
  } // end loop over elements
  // pack end-of-element-record
  pcbuff->packInt (LOADBALLANCER_END_DATA);
  
  return 1;
}


int
Domain::unpackMigratingData ( LoadBallancer* lb, ProcessCommunicator& pc) 
{
  // create temp space for dofManagers and elements
  // merging should be made by domain ?
  // maps of new dofmanagers and elements indexed by global number

  // we can put local dofManagers and elements into maps (should be done before unpacking)
  // int nproc=this->giveEngngModel()->giveNumberOfProcesses();
  int myrank=this->giveEngngModel()->giveRank();
  int iproc = pc.giveRank();
  int _mode, _globnum, _type;
  classType _etype;
  IntArray _partitions, local_partitions;
  //LoadBallancer::DofManMode dmode;
  DofManager* dofman;

  //  **************************************************
  //  Unpack migrating data to remote partition 
  //  **************************************************

  if (iproc == myrank) return 1; // skip local partition 
  // query process communicator to use
  ProcessCommunicatorBuff* pcbuff = pc.giveProcessCommunicatorBuff();
  ProcessCommDataStream pcDataStream (pcbuff);

  pcbuff->unpackInt (_type);
  // unpack dofman data
  while (_type != LOADBALLANCER_END_DATA) {
    _etype = (classType) _type;
    pcbuff->unpackInt (_mode);
    switch (_mode) {
    case LoadBallancer::DM_Remote: 
      // receiving new local dofManager
      pcbuff->unpackInt (_globnum);

      if (dmanMap.find(_globnum) != dmanMap.end()) { 
	// dofman is already available -> update only
        dofman = dmanMap[_globnum];
      } else { // data not available -> mode should be SharedUpdate
        dofman = CreateUsrDefDofManagerOfType (_etype, 0, this);
      }
      dofman->setGlobalNumber(_globnum);
      // unpack dofman state (this is the local dofman, not available on remote)
      dofman->restoreContext (&pcDataStream, CM_Definition | CM_State| CM_UnknownDictState);
      // unpack list of new partitions
      pcbuff->unpackIntArray (_partitions);
      dofman->setPartitionList (_partitions);
      dofman->setParallelMode (DofManager_local);
      dmanMap[_globnum] = dofman;
      break;

    case LoadBallancer::DM_Shared:      
      // receiving new shared dofManager, that was local on sending partition
      // should be received only once (from partition where was local)
      pcbuff->unpackInt (_globnum);

      if (dmanMap.find(_globnum) != dmanMap.end()) { 
	// dofman is already available -> update only
        dofman = dmanMap[_globnum];
      } else { // data not available -> mode should be SharedUpdate
        dofman = CreateUsrDefDofManagerOfType (_etype, 0, this);
      }
      dofman->setGlobalNumber(_globnum);
      // unpack dofman state (this is the local dofman, not available on remote)
      dofman->restoreContext (&pcDataStream, CM_Definition | CM_State| CM_UnknownDictState);
      // unpack list of new partitions
      pcbuff->unpackIntArray (_partitions);
      dofman->setPartitionList (_partitions);
      dofman->setParallelMode (DofManager_shared);
#if __VERBOSE_PARALLEL
      fprintf (stderr, "[%d] received Shared new dofman [%d]\n", myrank, _globnum);
#endif
      dmanMap[_globnum] = dofman;
      break;
       
    default:
      OOFEM_ERROR ("LoadBallancer::unpackMigratingData: unexpected dof manager type");
    }
    // get next type record
    pcbuff->unpackInt (_type);
    
  } ; // while (_type != LOADBALLANCER_END_DATA);
  
  // unpack element data
  Element* elem;
  do {
    pcbuff->unpackInt (_type);
    if (_type == LOADBALLANCER_END_DATA) break;
    _etype = (classType) _type;
    elem = CreateUsrDefElementOfType (_etype, 0, this);
    elem->restoreContext (&pcDataStream, CM_Definition | CM_State);
    recvElemList.push_back(elem);    
  } while (1);

  return 1;
}


// general algorithm common to all load ballancers
void 
Domain::migrateLoad (LoadBallancer* lb)
{
  int nproc=this->giveEngngModel()->giveNumberOfProcesses();
  int myrank=this->giveEngngModel()->giveRank();
  CommunicatorBuff cb (nproc, CBT_dynamic);
  Communicator com (this->giveEngngModel(), &cb, myrank, nproc, CommMode_Dynamic);

  // move existing dofmans and elements, that will be local on current partition,
  // into local map
  com.packAllData (this, lb, &Domain::packMigratingData);
  com.initExchange (MIGRATE_LOAD_TAG);
  
  // do something in between 
  this->initGlobalDofManMap ();
  this->deleteRemoteDofManagers (lb);
  this->deleteRemoteElements (lb);
  
  // receive remote data
  com.unpackAllData (this, lb, &Domain::unpackMigratingData);

  // Compress local data 

  AList<DofManager> *dofManagerList_new = new AList<DofManager>(0) ;
  AList<Element>    *elementList_new    = new AList<Element>(0) ;

  this->renumberDofManagers ();
  this->renumberDofManData () ;
  this->initializeNewDofManList (dofManagerList_new);

  this->compressElementData (elementList_new, lb);
  this->renumberElementData (lb);

  // 
  this->dofManagerList->clear(false); // not the data
  this->dofManagerList = dofManagerList_new;

  this->elementList->clear(false);
  this->elementList = elementList_new;

  // clean up
  this->dmanMap.clear();
  this->recvElemList.clear();

#if 1
  // debug print
  int i, j, nnodes=giveNumberOfDofManagers(), nelems=giveNumberOfElements();
  fprintf (stderr, "\n[%d] Nodal Table\n", myrank);
  for (i=1; i<=nnodes; i++) {
    if (giveDofManager(i)->giveParallelMode()==DofManager_local) 
      fprintf (stderr, "[%d]: %5d[%d] local\n", myrank, i, giveDofManager(i)->giveGlobalNumber());
    else if (giveDofManager(i)->giveParallelMode()==DofManager_shared) {
      fprintf (stderr, "[%d]: %5d[%d] shared ", myrank, i, giveDofManager(i)->giveGlobalNumber());
      for (j=1; j<=giveDofManager(i)->givePartitionList()->giveSize(); j++) {
        fprintf (stderr, "%d ", giveDofManager(i)->givePartitionList()->at(j));
      }
      fprintf (stderr, "\n");
    }
  }
  
  fprintf (stderr, "\n[%d] Element Table\n", myrank);
  for (i=1; i<=nelems; i++) {
    fprintf (stderr, "%5d {", i);
    for (j=1; j<=giveElement(i)->giveNumberOfDofManagers(); j++)
      fprintf (stderr, "%d ", giveElement(i)->giveDofManager(j)->giveNumber());
    fprintf (stderr, "}\n");
  }
#endif
  
}

void 
Domain::initGlobalDofManMap (bool forceinit)
{
  // initializes global dof man map according to domain dofman list

  if (forceinit || !dmanMapInitialized) {

    int key, idofman, ndofman = this->giveNumberOfDofManagers();
    DofManager* dofman;
    dmanMap.clear();
    
    for (idofman=1; idofman <= ndofman; idofman++) {
      dofman = this->giveDofManager (idofman);
      key = dofman->giveGlobalNumber();
      dmanMap[key] = dofman;
    }
  }
}


/* will delete those dofmanagers, that were sent to remote partition and are locally owned here
   so they are no longer necessary (those with state equal to DM_Remote and DM_SharedMerge)
   This will update domain DofManager list as well as global dmanMap and physically deletes the remote dofManager
*/
void
Domain::deleteRemoteDofManagers (LoadBallancer* lb)
{
  int i, ndofman =  this->giveNumberOfDofManagers();
  //LoadBallancer* lb = this->giveLoadBallancer();
  LoadBallancer::DofManMode dmode;
  DofManager* dman;
  int myrank=this->giveEngngModel()->giveRank();
  // loop over local nodes
  
  for (i = 1; i<= ndofman; i++) {
    dmode = lb->giveDofManState(i);
    if ((dmode == LoadBallancer::DM_Remote)) {
      // positive candidate found
      dmanMap.erase (this->giveDofManager (i)->giveGlobalNumber());
      //
      // this->deleteDofManager (i);  // delete and set entry to NULL
      //
      dman = dofManagerList->unlink (i);
      delete dman;
    } else if (dmode == LoadBallancer::DM_Shared) {
      dman = this->giveDofManager (i);
      dman->setPartitionList (*(lb->giveDofManPartitions(i)));
      dman->setParallelMode (DofManager_shared);
      if (!dman->givePartitionList()->findFirstIndexOf (myrank)) {
        dmanMap.erase (this->giveDofManager (i)->giveGlobalNumber());
        dman = dofManagerList->unlink (i);
        delete dman;
      }
    } else if (dmode == LoadBallancer::DM_Local) {
      IntArray _empty(0);
      dman = this->giveDofManager (i);
      dman->setPartitionList(_empty);
      dman->setParallelMode (DofManager_local);
    } else {
      OOFEM_ERROR ("Domain::deleteRemoteDofManagers: unknown dmode encountered");
    }
  }
}

/* will delete those elements, that were sent to remote partition and are locally owned here
   so they are no longer necessary (those with state equal to DM_Remote and DM_SharedMerge)
   This will update domain DofManager list as well as global dmanMap and physically deletes the remote dofManager
*/
void
Domain::deleteRemoteElements (LoadBallancer* lb)
{
  int i, nelem =  this->giveNumberOfElements();
  int myrank=this->giveEngngModel()->giveRank();
  //LoadBallancer* lb = this->giveLoadBallancer();
  Element* elem;

  // loop over local nodes
  
  for (i = 1; i<= nelem; i++) {
    if (lb->giveElementPartition(i) != myrank) {
      // positive candidate found
      // this->deleteElement (i);  // delete and set entry to NULL
      elem = elementList->unlink (i);
      delete (elem);
    }
  }
}



/*
  Assigns new local number (stored as dofmanager number, so it can be requested) 
  Assigns new local number to all dofManagers available in domanMap.
*/
void 
Domain::renumberDofManagers ()
{
  int _locnum;
  std::map<int, DofManager*>::iterator it;

  for (_locnum=1, it=dmanMap.begin(); it!=dmanMap.end(); it++) {
    it->second->setNumber(_locnum++);
  }
}


int
Domain::LB_giveUpdatedLocalNumber (int num, EntityRenumberingScheme scheme)
{
  if (scheme == ERS_DofManager) {
    DofManager* dm = this->giveDofManager (num);
    if (dm) {
      return dm->giveNumber();
    } else {
      _error2 ("LB_giveUpdatedLocalNumber: dofman %d moved to remote partition, updated number not available", num);
    }
  } else {
    _error ("LB_giveUpdatedLocalNumber: unsuported renumbering scheme");
  }
  return 0;
}

int
Domain::LB_giveUpdatedGlobalNumber (int num, EntityRenumberingScheme scheme)
{
  if (scheme == ERS_DofManager) {
    DofManager* dm = dmanMap[num];
    if (dm) {
      return dm->giveNumber();
    } else {
      _error2 ("LB_giveUpdatedGlobalNumber: dofman [%d] not available on local partition, updated number not available", num);
      return 0;
    }
  } else {
    _error ("LB_giveUpdatedGlobalNumber: unsuported renumbering scheme");
  }
  return 0;
}


void 
Domain::initializeNewDofManList (AList<DofManager>* dofManagerList)
{
  int _i, _size = dmanMap.size();
  std::map<int, DofManager*>::iterator it;
  dofManagerList->clear();
  dofManagerList->growTo(_size);

  for (_i=0, it=dmanMap.begin(); it!=dmanMap.end(); it++) {
    dofManagerList->put (++_i, it->second);
  }
}

void
Domain::compressElementData (AList<Element>* elementList, LoadBallancer* lb)
{
  int i, count, nelem = 0;
  int myrank=this->giveEngngModel()->giveRank();
  //LoadBallancer* lb = this->giveLoadBallancer();

  // determine real number of elements (excluding remote ones)
  
  for (i = 1; i<= this->giveNumberOfElements(); i++) {
    if ((lb->giveElementPartition(i) == myrank) && (this->giveElement(i) != NULL)) nelem++;
  }
  nelem+= recvElemList.size();
  
  elementList->clear(); elementList->growTo (nelem);
  count = 1;
  for (i = 1; i<= this->giveNumberOfElements(); i++) {
    if ((lb->giveElementPartition(i) == myrank) && (this->giveElement(i) != NULL)) 
      elementList->put(count++, this->giveElement(i));
  }
  std::list<Element*>::const_iterator it;
  for (it = recvElemList.begin(); it != recvElemList.end(); it++) {
    elementList->put(count++, *it);
  }
}

void 
Domain::renumberElementData (LoadBallancer* lb) {

  int i, nelems = this->giveNumberOfElements();
  int myrank = engineeringModel->giveRank();
  //LoadBallancer* lb = this->giveLoadBallancer();

  // loop first over local elements
  for (i = 1; i<= nelems; i++) {
    // skip remote ones
    if (lb->giveElementPartition(i) != myrank) continue;
    //  ask elelement to update numbering
    this->giveElement(i)->updateLocalNumbering (this, &Domain::LB_giveUpdatedLocalNumber); // should call LoadBallncer for translation
  }
  // now loop over received elemens (they have node numbers in global numbering)
  std::list<Element*>::const_iterator it;
  for (it = recvElemList.begin(); it != recvElemList.end(); it++) {
    (*it)->updateLocalNumbering (this, &Domain::LB_giveUpdatedGlobalNumber); // g_to_l
  }
}

/* renumber here the master node number for rigid and hanging dofs, etc;
   existing local nodes need mapping from old_local to new numbering,
   but received nodes need mapping from global to new numbering

   -> we need to keep the list of received nodes! (now they are only introduced into globally indexed dmanMap!);
*/
void
Domain::renumberDofManData () {
  int _i;
  std::map<int, DofManager*>::iterator it;

  for (_i=0, it=dmanMap.begin(); it!=dmanMap.end(); it++) {
    if (it->second->giveNumber() == 0) {
      // received dof manager -> we map global numbers to new local number
      it->second->updateLocalNumbering (this, &Domain::LB_giveUpdatedGlobalNumber); // g_to_l
    } else {
      // existing dof manager -> we map old local number to new local number
      it->second->updateLocalNumbering (this, &Domain::LB_giveUpdatedLocalNumber); // l_to_l
    }
  }
}

int
Domain::dofmanGlobal2local (int _globnum)
{
  if (dmanMap.find(_globnum) != dmanMap.end()) { 
    // dofman is already available -> update only
    return (dmanMap[_globnum]->giveNumber());
  } else return 0;
}

#endif
