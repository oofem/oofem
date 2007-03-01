/* $Header: /home/cvs/bp/oofem/oofemlib/src/compcol.C,v 1.5.4.1 2004/04/05 15:19:43 bp Exp $ */
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

#include "dss.h"
#include "error.h"

#ifdef __DSS_MODULE

#include "flotarry.h"
#include "engngm.h"
#include "domain.h"
#include "DSSolver.h"
#ifndef __MAKEDEPEND
#include <set>
#endif

DSSMatrix :: DSSMatrix (dssType _t): SparseMtrx()
{
  eDSSolverType _st = eDSSFactorizationLDLT;
  _dss = new DSSolver();
  _type = _t;
  if (_t == sym_LDL) _st=eDSSFactorizationLDLT;
  else if (_t==sym_LL) _st = eDSSFactorizationLLT;
  else if (_t==unsym_LU) _st = eDSSFactorizationLU;
  else OOFEM_ERROR ("DSSMatrix::DSSMatrix() -- unknown dssType"); 
  _dss->Initialize (0, _st);
  isFactorized = FALSE;
}


DSSMatrix::DSSMatrix (dssType _t, int n) : SparseMtrx(n,n)
{
  eDSSolverType _st = eDSSFactorizationLDLT;
  _dss = new DSSolver();
  _type = _t;
  if (_t == sym_LDL) _st=eDSSFactorizationLDLT;
  else if (_t==sym_LL) _st = eDSSFactorizationLLT;
  else if (_t==unsym_LU) _st = eDSSFactorizationLU;
  else OOFEM_ERROR ("DSSMatrix::DSSMatrix() -- unknown dssType"); 
  _dss->Initialize (0, _st);
  isFactorized = FALSE;
}   

DSSMatrix::~DSSMatrix ()
{
  delete _dss;
}

/*****************************/
/*  Copy constructor         */
/*****************************/

DSSMatrix::DSSMatrix(const DSSMatrix &S)
{
  OOFEM_ERROR ("DSSMatrix::DSSMatrix(const DSSMatrix &S) -- not implemented"); 
}



/***************************/
/* Assignment operator...  */
/***************************/
/*
DSSMatrix& DSSMatrix::operator=(const DSSMatrix &C)  
{
  OOFEM_ERROR ("DSSMatrix::operator= -- not implemented"); 
  return *this;
}
*/


SparseMtrx* DSSMatrix::GiveCopy () const
{
  OOFEM_ERROR ("DSSMatrix::GiveCopy -- not implemented"); 
  return NULL;
}


void DSSMatrix::times (const FloatArray& x, FloatArray& answer) const
{
  OOFEM_ERROR ("DSSMatrix::times -- not implemented"); 
  return ;
}

void DSSMatrix::times (double x) 
{
  OOFEM_ERROR ("DSSMatrix::times -- not implemented"); 
  return ;
}

int DSSMatrix::buildInternalStructure (EngngModel*eModel, int di, EquationID ut) 
{
 IntArray  loc;
 Domain* domain = eModel->giveDomain(di);
 int neq = eModel -> giveNumberOfDomainEquations (di,ut);
 int nelem = domain -> giveNumberOfElements() ;
 int i,ii,j,jj,n, indx;
 Element* elem;
 // allocation map 
 std::vector< std::set<int> > columns(neq);

 int nz_ = 0;

 for (n=1 ; n<=nelem ; n++) {  
  elem = domain -> giveElement(n);
  elem -> giveLocationArray (loc, ut) ;

  for (i=1 ; i <= loc.giveSize() ; i++) {
   if ((ii = loc.at(i))) {
    for (j=1; j <= loc.giveSize() ; j++) {
     if ((jj=loc.at(j))) {
      columns[jj-1].insert(ii-1);
     }
    }
   }
  }
 }

 for (i=0; i<neq; i++) nz_ += columns[i].size();
  
 unsigned long rowind_[nz_];
 unsigned long colptr_[neq+1];
 indx = 0;

 std::set<int>::iterator pos;
 for (j=0; j<neq; j++) { // column loop
  colptr_[j] = indx;
  for (pos = columns[j].begin(); pos != columns[j].end(); ++pos) { // row loop
   rowind_[indx++] = *pos;
  }
 }
 colptr_[neq] = indx;

 SparseMatrixF sm (neq, NULL, rowind_, colptr_, 0,0,true);
 int bsize = eModel->giveDomain(1)->giveDefaultNodeDofIDArry().giveSize();
 _dss -> SetMatrixPattern (&sm, bsize) ;
 _dss -> PreFactorize ();
    
 OOFEM_LOG_DEBUG ("DSSMatrix info: neq is %d, bsize is %d\n",neq,nz_);

 // increment version
 this->version++;
 
 return TRUE;
}


int DSSMatrix::assemble (const IntArray& loc, const FloatMatrix& mat)
{
 int        i,j,ii,jj,dim ;
 
#  ifdef DEBUG
 dim = mat.giveNumberOfRows() ;
 if (dim != loc.giveSize()) 
   OOFEM_ERROR ("CompCol::assemble : dimension of 'k' and 'loc' mismatch") ;
#  endif
 
 dim = mat.giveNumberOfRows() ;

 if (_type == unsym_LU) {
   for (j=1; j<=dim; j++) {
     jj = loc.at(j);
     if(jj){
       for (i=1; i<=dim; i++) {
         ii = loc.at(i);
         if(ii){
           _dss->ElementAt(ii-1,jj-1) += mat.at(i,j); 
         }
       }
     }
   }
 } else { // symmetric pattern
   for (j=1; j<=dim; j++) {
     jj = loc.at(j);
     if(jj){
       for (i=j; i<=dim; i++) {
         ii = loc.at(i);
         if(ii){
           _dss->ElementAt(ii-1,jj-1) += mat.at(i,j); 
         }
       }
     }
   }
 }
 
 // increment version
 this->version++;
 
 return 1;

}

int DSSMatrix::assemble (const IntArray& rloc, const IntArray& cloc, const FloatMatrix& mat)
{

 int        i,j,ii,jj,dim1,dim2 ;

 // this->checkSizeTowards(rloc, cloc);
  
 dim1 = mat.giveNumberOfRows() ;
 dim2 = mat.giveNumberOfColumns() ;
 for (i=1 ; i<= dim1; i++) { 
  ii = rloc.at(i);
  if (ii) for (j=1 ; j<= dim2; j++) {
   jj = cloc.at(j);
   if (jj) _dss->ElementAt(ii-1,jj-1) += mat.at(i,j);
  }
 }
 
 // increment version
 this->version++;
 
 return 1;
 
}

SparseMtrx* DSSMatrix::zero ()
{
  _dss->LoadZeros();
  
  // increment version
  this->version++;
  isFactorized = FALSE ;
  
  return this;
}

SparseMtrx* DSSMatrix::factorized ()
{
 if (isFactorized) return this;

 _dss->ReFactorize();
 isFactorized = TRUE ;
 return this ;
}

void
DSSMatrix::solve (FloatArray* b, FloatArray* x)
{
  x->resize(b->giveSize());
  _dss->Solve(x->givePointer(), b->givePointer());
}

/*
void DSSMatrix::toFloatMatrix (FloatMatrix& answer) const
{}
*/
/*
void DSSMatrix::printYourself () const
{}
*/
/*********************/
/*   Array access    */
/*********************/

double& DSSMatrix::at (int i, int j) 
{
  
 // increment version
 this->version++;
 return _dss->ElementAt(i-1,j-1);
}


double  DSSMatrix::at (int i, int j) const
{
 return _dss->ElementAt(i-1,j-1);
}

double DSSMatrix::operator()(int i, int j)  const
{
 return _dss->ElementAt(i,j);
}

double& DSSMatrix::operator()(int i, int j) 
{        
 // increment version
 this->version++;
 return _dss->ElementAt(i,j);
}

#else // ifndef __DSS_MODULE

double DSS__zero ;

DSSMatrix::DSSMatrix (dssType _t,int n) : SparseMtrx() {
  OOFEM_ERROR ("DSSMatrix: can't create, DSS support not compiled");
}

DSSMatrix::DSSMatrix (dssType _t) : SparseMtrx() {
  OOFEM_ERROR ("DSSMatrix: can't create, DSS support not compiled");
}

DSSMatrix::DSSMatrix(const DSSMatrix &S) {
  OOFEM_ERROR ("DSSMatrix: can't create, DSS support not compiled");
}

DSSMatrix::~DSSMatrix () {}

SparseMtrx* DSSMatrix::GiveCopy () const {
  OOFEM_ERROR ("DSSMatrix: can't create, DSS support not compiled");
}

void DSSMatrix::times (const FloatArray& x, FloatArray& answer) const {}
void DSSMatrix::times (double x) {}
int DSSMatrix::buildInternalStructure (EngngModel*, int, EquationID) {return 0;}
int DSSMatrix::assemble (const IntArray& loc, const FloatMatrix& mat) {return 0;}
int DSSMatrix::assemble (const IntArray& rloc, const IntArray& cloc, const FloatMatrix& mat) {return 0;}
SparseMtrx* DSSMatrix::factorized () {return  NULL;}
void DSSMatrix::solve (FloatArray* b, FloatArray* x) {}
SparseMtrx* DSSMatrix::zero () {return NULL;}
double& DSSMatrix::at (int i, int j) {return DSS__zero;}
double DSSMatrix::at (int i, int j) const {return 0.;}
double DSSMatrix::operator() (int i, int j) const {return 0.;};               
double& DSSMatrix::operator() (int i, int j) {return DSS__zero;}               

#endif       


