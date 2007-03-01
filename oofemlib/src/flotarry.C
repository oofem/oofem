/* $Header: /home/cvs/bp/oofem/oofemlib/src/flotarry.C,v 1.17.4.1 2004/04/05 15:19:43 bp Exp $ */
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
/*
 The original idea for this class comes from 
  Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 PhD Thesis, EPFL, Lausanne, 1992.
*/

//   file FLOTARRY.CC

#include "flotarry.h"
#include "intarray.h"
#include "flotmtrx.h"
#include "mathfem.h"
#include "error.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <string.h>
#endif

#ifdef __PARALLEL_MODE
#include "combuff.h"
#endif

FloatArray :: FloatArray (int n)
   // Constructor : creates an array of size n (filled with garbage).
{
   allocatedSize = size = n ;
   if (size)
      values = allocDouble(size) ;
   else
      values = NULL ;
}

FloatArray :: FloatArray (const FloatArray & src) 
{
// copy constructor

 double * srcVal;

 allocatedSize = size = src.size ;
 if (size)
  values = allocDouble(size) ;
 else
  values = NULL ;
 
 srcVal = src.givePointer();
 for (int i=0; i< size; i++) this->values[i]=srcVal[i];
}

FloatArray&
FloatArray :: operator= (const FloatArray& src)
{
 // assignment: cleanup and copy
 double *srcVal;

 if (this != &src) { // beware of s=s;
  this->resize (src.size);
  
  srcVal = src.givePointer();
  for (int i=0; i< size; i++) this->values[i]=srcVal[i];
 }
 return *this;
}


void
FloatArray :: add (const FloatArray& b)
   // Performs the operation a=a+b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (b.giveSize()==0)
      return  ;

   if (! size)                                // null-sized array
   {resize (b.size); zero();}

  if (size != b.size) {                  // unmatching sizes
   OOFEM_ERROR3 ("FloatArray :: add :  dimension mismatch in a[%d]->add(b[%d])\n", size,b.size) ;
  }

   p1 = values ;
   p2 = b.values ;
   i  = size ;
   while (i--)
   *p1++ += *p2++ ;
   return  ;
}


void
FloatArray :: add (const double factor, const FloatArray& b)
   // Performs the operation a=a+factor*b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (b.giveSize()==0)
      return  ;

   if (! size)                                // null-sized array
		 {resize (b.size); zero();}

	 if (size != b.size) {                  // unmatching sizes
		 OOFEM_ERROR3 ("FloatArray :: add :  dimension mismatch in a[%d]->add(b[%d])\n", size,b.size) ;
   }

   p1 = values ;
   p2 = b.values ;
   i  = size ;
   while (i--)
		 *p1++ += factor * (*p2++) ;
   return  ;
}

void  FloatArray :: substract (const FloatArray& src)
   // Performs the operation a=a-src, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of src. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (src.giveSize()==0)
      return  ;

   if (! size) {                              // null-sized array
   this->resize (size = src.size) ;
   this->zero();
  }

#  ifdef DEBUG
  if (size != src.size) {                  // unmatching sizes
    OOFEM_ERROR3 ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n", size,src.size) ;
  }
#  endif

   p1 = values ;
   p2 = src.values ;
   i  = size ;
   while (i--)
      *p1++ -= *p2++ ;
   return  ;
}



double  dotProduct (const FloatArray& p1, const FloatArray& p2, register int i)
   // A non-member function. Returns the dot product of the first 'i' coef-
   // ficienst of the two arrays P1 and P2. This method applies to many
   // situations, eg row*column products with matrices.
{
   double answer ;

   answer = 0. ;
   while (i) {
   answer += p1.at(i) * p2.at(i) ;
   i--;
  }
   return answer ;
}

void
FloatArray :: beSubArrayOf (const FloatArray& src, const IntArray &indx)
//
// returns subVector from receiver
// subVector size will be indx max val size
// and on i-th position of subVector will be this->at(indx->at(i))
//
{
 //FloatArray *answer;
 int i,ii,n,isize;

 n = indx.giveSize();

 if (src.size != n) {OOFEM_ERROR ("FloatArray :: beSubArrayOf - size mismatch");}

 for (isize=0,i=1; i<=n;i++) if(indx.at(i)>isize) isize = indx.at(i);
 //answer = new FloatArray (isize);
 this->resize (isize);
 for (i=1; i<= n; i++) {
  ii = indx.at(i);
  if (ii > 0) this->at(ii) = src.at(i);
 }
 return ;
}


void 
FloatArray :: addSubVector (const FloatArray& src, int si)
{
 int i, reqSize, n = src.giveSize();

 si --;
 reqSize = si+n;
 if (this->giveSize() < reqSize) this->resize(reqSize);

 for (i=1; i<= n; i++) this->at(si+i) += src.at(i);
}



void
FloatArray :: beVectorProductOf (const FloatArray& v1, const FloatArray& v2)
//
// computes vector product v1 x v2
// and stores result into receiver
//
{
 // check proper bunds
 if ((v1.giveSize() != 3) || (v2.giveSize() != 3)) 
  {
   OOFEM_ERROR (" FloatArray::VectorProduct : size mismatch, size is not equal to 3") ;
  }
 this->resize (3);
 
 this->at(1) = v1.at(2)*v2.at(3)-v1.at(3)*v2.at(2);
 this->at(2) = v1.at(3)*v2.at(1)-v1.at(1)*v2.at(3);
 this->at(3) = v1.at(1)*v2.at(2)-v1.at(2)*v2.at(1);
 
 return ;
}

double
FloatArray :: distance (const FloatArray &from) const
//
// returns distance between receiver and from from
// computed using generalized pythagora formulae
//
{
  double dist = 0.;
  
/*
  if (size != from.giveSize())  {
    printf ("error in FloatArray::distance size mismatch (%d) is not equal to (%d) \n",
            size,from.giveSize()) ;
    exit(0) ;}
*/
  int s = min (size, from.giveSize());
  for (int i = 1; i<= s; i++)
    dist += (this->at(i)-from.at(i))*(this->at(i)-from.at(i));
  
  return sqrt(dist);
 
} 


/*******************************************************************************/

FloatArray*  FloatArray :: add (FloatArray* b)
   // Performs the operation a=a+b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (!b || b->giveSize()==0)
      return this ;

   if (! size) {                              // null-sized array
      size   = b -> size ;
      values = allocDouble(size) ;}

#  ifdef DEBUG
      if (size != b->size) {                  // unmatching sizes
        OOFEM_ERROR3 ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n", size,b->size) ;
      }
#  endif

   p1 = values ;
   p2 = b -> values ;
   i  = size ;
   while (i--)
      *p1++ += *p2++ ;
   return this ;
}


FloatArray*  FloatArray :: substract (FloatArray* b)
   // Performs the operation a=a-b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (!b || b->giveSize()==0)
      return this ;

   if (! size) {                              // null-sized array
      size   = b -> size ;
      values = allocDouble(size) ;}

#  ifdef DEBUG
      if (size != b->size) {                  // unmatching sizes
        OOFEM_ERROR3 ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n", size,b->size) ;
      }
#  endif

   p1 = values ;
   p2 = b -> values ;
   i  = size ;
   while (i--)
      *p1++ -= *p2++ ;
   return this ;
}



void  FloatArray :: assemble (const FloatArray& fe, const IntArray& loc)
   // Assembles the array fe (typically, the load vector of a finite
   // element) to the receiver, using loc as location array.
{
   int i,ii,n ;

#  ifdef DEBUG
      if ((n=fe.giveSize()) != loc.giveSize()) {
        OOFEM_ERROR ("FloatArray::assemble : dimensions of 'fe' and 'loc' mismatch") ;
      }
      this -> checkSizeTowards(loc) ;
#  endif

   n = fe.giveSize() ;
   for (i=1 ; i<=n ; i++) {
      ii = loc.at(i) ;
      if (ii)                                  // if non 0 coefficient,
  this->at(ii) += fe.at(i) ; }         // then assemble
}


#ifdef DEBUG
double&  FloatArray :: at (int i)
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}

double  FloatArray :: at (int i) const
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}
#endif


#ifdef DEBUG
void  FloatArray :: checkBounds(int i) const
   // Checks that the receiver's size is not smaller than 'i'.
{
  if (i<=0) 
    OOFEM_ERROR2 ("FloatArray :: checkBounds : array error on index : %d <= 0 \n",i) ;
  if (i>size) 
    OOFEM_ERROR3 ("FloatArray :: checkBounds : array error on index : %d > %d \n",i,size) ;
}
#endif

FloatArray*
FloatArray :: GiveSubArray (IntArray *indx)
//
// returns subVector from receiver
// subVector size will be indx max val size
// and on i-th position of subVector will be this->at(indx->at(i))
//
{
 FloatArray *answer;
 int i,ii,n,isize;

 n = indx->giveSize();

 if (size != n) OOFEM_ERROR ("FloatArray :: GiveSubArray : size mismatch");

 for (isize=0,i=1; i<=n;i++) if(indx->at(i)>isize) isize = indx->at(i);
 answer = new FloatArray (isize);
 for (i=1; i<= n; i++) {
  ii = indx->at(i);
  if (ii > 0)answer->at(ii) = this->at(i);
 }
 return answer;
}

void  FloatArray :: checkSizeTowards (const IntArray& loc)
   // Expands the receiver if loc points to coefficients beyond the size of
   // the receiver.
{
   int i,n,high ;

   high = 0 ;
   n    = loc.giveSize() ;
   for (i=1 ; i<=n ; i++)
      high = max (high,(loc.at(i))) ;
   if (high > size)                             // receiver must be expanded
      this -> resize (high) ;
}


void  FloatArray :: resize (int n, int allocChunk)
   // Expands the receiver up to size n (n is assumed larger than 'size').
   // Initializes all new coefficients to zero.
{
   register int i ;
   double       *newValues,*p1,*p2 ;

   if (n <= allocatedSize) {size = n; return;}

  if (allocChunk < 0) allocChunk = 0;
   newValues = allocDouble(n + allocChunk) ;

   p1 = values ;
   p2 = newValues ;
   i  = size ;
   while (i--)
      *p2++ = *p1++ ;

   if (values)
   freeDouble (values) ;
   values = newValues ;
   allocatedSize = n + allocChunk;
  size   = n ;
 }


void  FloatArray :: hardResize (int n)
   // Realocates the receiver with new size.
   // Initializes all new coefficients to zero.
{
   register int i ;
   double       *newValues,*p1,*p2 ;

   // if (n <= allocatedSize) {size = n; return;}

   newValues = allocDouble(n) ;

   p1 = values ;
   p2 = newValues ;
   i  = min(size, n) ;
   while (i--)
      *p2++ = *p1++ ;

   if (values)
      freeDouble (values) ;
   values = newValues ;
   allocatedSize = size   = n ;
}



int  FloatArray :: containsOnlyZeroes () const
   // Returns True if all coefficients of the receiver are 0, else returns
   // False.
{
   register int i ;
   double       *p ;

   p = values ;
   i = size ;
   while (i--)
      if (*p++ != 0.)
  return FALSE ;

   return TRUE ;
}



void FloatArray :: zero ()
// zeroes all values to zero
{
  int i;
  if (values)
  for (i=0; i< size; i++) values[i]=0.;
  // return *this ;
}


void  FloatArray :: beProductOf (const FloatMatrix& aMatrix, const FloatArray& anArray)
   // Stores the product of aMatrix * anArray in to receiver
{
  int         i,j, nColumns, nRows ;
  double      sum ;
  //FloatArray* answer ;
  
#  ifdef DEBUG
  if ((nColumns=aMatrix.giveNumberOfColumns()) - anArray.giveSize()) {
    OOFEM_ERROR ("FloatArray :: beProductOf : dimension mismatch");
  }
#  endif
  
  nColumns = aMatrix.giveNumberOfColumns();
  this->resize (nRows = aMatrix.giveNumberOfRows()) ;
  for (i=1 ; i<=nRows ; i++) {
      sum = 0. ;
      for (j=1 ; j<=nColumns ; j++)
    sum += aMatrix.at(i,j) * anArray.at(j) ;
      this->at(i) = sum ;}
   return ;
}

void  FloatArray :: beTProductOf (const FloatMatrix& aMatrix, const FloatArray& anArray)
   // Stores the product of aMatrix^T * anArray in to receiver
{
   int         i,j, nColumns, nRows ;
   double      sum ;
   //FloatArray* answer ;

#  ifdef DEBUG
   if ((nColumns=aMatrix.giveNumberOfRows()) - anArray.giveSize()) {
     OOFEM_ERROR ("FloatArray :: beTProductOf : dimension mismatch");
   }
#  endif

  nColumns = aMatrix.giveNumberOfRows();
   this->resize (nRows = aMatrix.giveNumberOfColumns()) ;
   for (i=1 ; i<=nRows ; i++) {
      sum = 0. ;
      for (j=1 ; j<=nColumns ; j++)
    sum += aMatrix.at(j,i) * anArray.at(j) ;
      this->at(i) = sum ;}
   return ;
}


FloatArray*  FloatArray :: negated ()
   // Switches the sign of every coefficient of the receiver. Returns the
   // receiver.
{
   register int i ;
   double       x ;
   double*      p ;

   i = size ;
   p = values ;
   while (i--) {
      x    = - *p ;
      *p++ = x ;}
   return this ;
}


void  FloatArray :: printYourself () const
   // Prints the receiver on screen.
{
   printf ("FloatArray of size : %d \n",size) ;
   for (int i=1 ; i<=size ; ++i)
      printf ("%10.3e  ",this->at(i)) ;
   printf ("\n") ;
}

FloatArray*  FloatArray :: setValuesToZero ()
// zeroes all values to zero
{
  int i;
  if (values)
  for (i=0; i< size; i++) values[i]=0.;
  return this ;
}


FloatArray*  FloatArray :: beCopyOf (FloatArray*arry)
// Returns the receiver initialized according to array arry
// if arry is NULL nothing done
{
  int i;
  double *toVal;
  if (arry) {
  if (size != arry->giveSize()) {
  freeDouble (values);
    size = arry->giveSize();
  values = allocDouble(size) ;
  }
  
  toVal = arry->givePointer();
  for (i=0; i< size; i++) values[i]=toVal[i];
  return this ;
  } else return this;
}
  


void
FloatArray :: rotatedWith (FloatMatrix& r, char mode)
   // Returns the receiver 'a' rotated according the change-of-base matrix r.
   // If mode = 't', the method performs the operation  a = r(transp) * a .
   // If mode = 'n', the method performs the operation  a = r * a .
{
   FloatMatrix  rot ;
   FloatArray   rta ;

   if (mode == 't') {
      rot.beTranspositionOf(r) ;
   rta.beProductOf (rot, *this);
  } else if (mode == 'n') {
   rta.beProductOf (r, *this);
  } else {
    OOFEM_ERROR ("FloatArray :: rotatedWith: unsupported mode");
  }

  *this = rta;
}


FloatArray*  FloatArray :: rotatedWith (FloatMatrix* r, char mode)
   // Returns the receiver 'a' rotated according the change-of-base matrix r.
   // If mode = 't', the method performs the operation  a = r(transp) * a .
   // If mode = 'n', the method performs the operation  a = r * a .
{
   // double       *p1,*p2 ;
   FloatMatrix  *rot = NULL ;
   FloatArray   *rta ;

   if (mode == 't')
   rot = r -> GiveTransposition() ;
   else if (mode == 'n') 
   rot = r ;
  else {
   OOFEM_ERROR ("FloatArray :: rotatedWith: unsupported mode");
  }

   rta = rot -> Times(this) ;

/*
   p1 = values ;
   p2 = rta -> values ;

   i  = size ;
   while (i--)
      *p1++ = *p2++ ;
*/

  *this = *rta;

   if (mode == 't')
      delete rot ;
   delete rta ;
   return this ;
}


FloatArray*  FloatArray :: times (double factor)
   // Multiplies every coefficient of the receiver by factor. Answers the
   // modified receiver.
{
   register int i ;
   double*      p ;

   p = values ;
   i = size ;
   while (i--)
      *(p++) *= factor ;
   return this ;
}


FloatArray*  FloatArray :: Times (double factor) const
   // Returns a new array, whose components are those of the receicer, times
   // factor.
{
   register int i ;
   double       *p1,*p2 ;
   FloatArray*  answer ;

   answer = new FloatArray(size) ;
   p1     = values ;
   p2     = answer -> values ;
   i      = size ;
   while (i--)
      *p2++ = factor * (*p1++) ;
   return answer ;
}


double  dotProduct (double* P1, double* P2, register int i)
   // A non-member function. Returns the dot product of the first 'i' coef-
   // ficienst of the two arrays P1 and P2. This method applies to many
   // situations, eg row*column products with matrices.
{
   double answer ;

   answer = 0. ;
   while (i--)
      answer += *P1++ * *P2++ ;
   return answer ;
}


double
FloatArray :: distance (FloatArray *from)
//
// returns distance between receiver and from from
// computed using generalized pythagora formulae
//
{
 double dist = 0.;

 if (size != from-> giveSize())  {
   OOFEM_ERROR3 ("FloatArray::distance :  size mismatch (%d) is not equal to (%d)", size,from->giveSize()) ;
 }
 
 for (int i = 1; i<= size; i++)
  dist += (this->at(i)-from->at(i)) * (this->at(i)-from->at(i));
 
 return sqrt(dist);
 
} 

FloatArray*
FloatArray :: VectorProduct (FloatArray* v2)
//
// computes vector product this x v2
// and return newly allocated result
//
{
 // check proper bunds
 if ((size != 3) || (v2->giveSize() != 3)) 
  {
    OOFEM_ERROR ("FloatArray::VectorProduct : size mismatch, size is not equal to 3") ;
  }
 FloatArray *answer = new FloatArray(3);
 
 answer->at(1) = this->at(2)*v2->at(3)-this->at(3)*v2->at(2);
 answer->at(2) = this->at(3)*v2->at(1)-this->at(1)*v2->at(3);
 answer->at(3) = this->at(1)*v2->at(2)-this->at(2)*v2->at(1);

 return answer;
}


FloatArray*
FloatArray :: normalize ()
//
// normalizes receiver to have norm equal to 1.0
//
{
  double norm = this->computeNorm();
  if (norm < 1.e-80) {
    OOFEM_ERROR ("FloatArray::normalize : cannot norm receiver, norm is too small");
  }
  
  this->times(1./norm);
  
 return this;
}


double
FloatArray :: computeNorm()
{
 int i;
 double norm = 0.;
 for (i=1; i<= size; i++)
  norm += this->at(i)*this->at(i);
 
 return sqrt (norm);
}
 


contextIOResultType  FloatArray :: storeYourself (FILE* stream)
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 succes
//              =0 file i/o error
{
  int type_id = FloatArrayClass;
  // write class header
  if (fwrite(&type_id,sizeof(int),1,stream)!= 1) return CIO_IOERR;
  // write size 
  if (fwrite(&size,sizeof(int),1,stream)!= 1) return CIO_IOERR;
  // write raw data
 if (size) {
  if (fwrite(values,sizeof(double),size,stream) != (size_t) size) return CIO_IOERR;
 }
  // return result back
  return CIO_OK;
}

contextIOResultType  FloatArray :: restoreYourself (FILE* stream)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id od class id is not correct
{
  int class_id;
  // read class header
  if (fread(&class_id,sizeof(int),1,stream) != 1) return CIO_IOERR;
  if (class_id != FloatArrayClass) return CIO_BADVERSION;
  // read size 
  if (fread(&size,sizeof(int),1,stream) != 1) return CIO_IOERR;
  if (values!=NULL) freeDouble(values);
  if (size) {
    values = allocDouble(size) ;
  allocatedSize = size;
  // read raw data
  if (fread(values,sizeof(double),size,stream) != (size_t)size) return CIO_IOERR;
  } else {
    values = NULL ;
  allocatedSize = 0;
 }

  // return result back
  return CIO_OK;
}


#ifdef __PARALLEL_MODE
int 
FloatArray :: packToCommBuffer (CommunicationBuffer& buff) const
{
 int result = 1;
 // pack size
 result &= buff.packInt (size);
 // pack data
 result &= buff.packDoubleArray (this->values, size);
 
 return result;
}

int 
FloatArray :: unpackFromCommBuffer (CommunicationBuffer& buff)
{
 int newSize, result = 1;
 // unpack size
 result &= buff.unpackInt (newSize);
 // resize yourself
 this->resize(newSize);
 result &= buff.unpackDoubleArray (this->values, newSize);

 return result;
}

int 
FloatArray :: givePackSize (CommunicationBuffer& buff) const
{
 return buff.giveDoubleVecPackSize (1) + buff.giveDoubleVecPackSize (this->size);
}

#endif


#ifdef IML_COMPAT

FloatArray & 
FloatArray::operator=(const double& val)
{
 for (int i=0; i< size; i++) values[i] = val;
 return *this;
}

FloatArray& operator*=(FloatArray &x, const double &a)
{
 int N = x.giveSize();
 for (int i=0;i<N;i++)
  x(i) *= a;
 return x;
}


FloatArray operator*(const double &a, const FloatArray &x)
{
 int N = x.giveSize();
 FloatArray result(N);
 for (int i=0;i<N;i++)
  result(i) = x(i)*a;
 return result;
}

FloatArray operator*(const FloatArray &x, const double &a)
{
 // This is the other commutative case of vector*scalar.
 // It should be just defined to be
 // "return operator*(a,x);"
 // but some compilers (e.g. Turbo C++ v.3.0) have trouble
 // determining the proper template match.  For the moment,
 // we'll just duplicate the code in the scalar * vector 
 // case above.

 int N = x.giveSize();
 FloatArray result(N);
 for (int i=0;i<N;i++)
  result(i) = x(i)*a;
 return result;
 
}

FloatArray operator+(const FloatArray &x, const FloatArray &y)
{
 int N = x.giveSize();
 if (N != y.giveSize())
  {
   OOFEM_ERROR ("loatArray operator+ : incompatible vector lengths");
  }
 
 FloatArray result(N);
 for (int i=0;i<N; i++)
  result(i) = x(i) + y(i);
 return result;
}
          
FloatArray operator-(const FloatArray &x, const FloatArray &y)
{
 int N = x.giveSize();
 if (N != y.giveSize())
  {
   OOFEM_ERROR ("FloatArray operator- : incompatible vector lengths");
  }
 
 FloatArray result(N);
 for (int i=0;i<N; i++)
  result(i) = x(i) - y(i);
 return result;
}
          

FloatArray& operator+=(FloatArray &x, const FloatArray &y)
{
 int N = x.giveSize();
 if (N != y.giveSize())
  {
    OOFEM_ERROR ("FloatArray& operator+= : incompatible vector lengths");
  }
      
 for (int i=0;i<N; i++)
  x(i) += y(i);
 return x;
}
          
      
FloatArray& operator-=(FloatArray &x, const FloatArray &y)
{
 int N = x.giveSize();
 if (N != y.giveSize())
  {
   OOFEM_ERROR ("FloatArray& operator-= : incompatible vector lengths");
  }
 
 for (int i=0;i<N; i++)
  x(i) -= y(i);
 return x;
}
          

double dot(const FloatArray &x, const FloatArray &y)
{
        
  //  Check for compatible dimensions:
  if (x.giveSize() != y.giveSize())
  {
   OOFEM_ERROR ("dot : incompatible dimensions") ;
  }
 
 double temp =  0;
 for (int i=0; i<x.giveSize();i++)
  temp += x(i)*y(i);
 return temp;
}

double norm(const FloatArray &x)
{
      double temp = dot(x,x);
      return sqrt(temp);
}

#endif


