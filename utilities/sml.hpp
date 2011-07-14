/* Copyright (c) 2009, NICTA 
 * All rights reserved. 
 * 
 * The contents of this file are subject to the Mozilla Public License 
 * Version 1.1 (the "License"); you may not use this file except in 
 * compliance with the License. You may obtain a copy of the License at 
 * http://www.mozilla.org/MPL/ 
 * 
 * Software distributed under the License is distributed on an "AS IS" 
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the 
 * License for the specific language governing rights and limitations 
 * under the License. 
 * 
 * Authors: Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *          S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
 *
 * Created: (02/11/2007) 
 *
 * Last Updated: (14/11/2007)   
 */

#ifndef _SML_HPP_
#define _SML_HPP_

#include "bmrmexception.hpp"
#include <iostream>
#include <cstring>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/storage.hpp>


namespace ublas = boost::numeric::ublas;


namespace SML {

  const double ZERO_EPS = 1e-20;  
  const double INFTY = 1e30;    // single precision: ~\pm 10^{38.53}; double precision: ~\pm 10^{308.25}
  
  enum VECMAT {VECTOR, MATRIX};
  enum COLROW {COLUMN,ROW};
  enum DENSPA {DENSE,SPARSE};
  enum TYPE   {DENSE_VECTOR,SPARSE_VECTOR,DENSE_MATRIX,SPARSE_MATRIX,SPARSE_MATRIX_ROW,DENSE_MATRIX_ROW};     // combination of dense/sparse and vector/matrix
  
   // absolute value function
   template<typename T>
   inline T abs(T x){ return (x>0) ? (x) : (-x); }
   
   // sign function
   template<typename T>
   inline int sgn(T x){ return (x>0) ? 1 : ((x<0) ? -1 : 0); }
   
   // square
   template<typename T>
   inline T sqr(T x){ return x*x;}

} // namespace SML



/**
 * Simple Matrix Library (SML)
 *
 * This is the dark underbelly of the solver framework. The aim of this
 * class is to opaquely perform matrix operations without the calling
 * ever needing to know what the representation of underlying matrices.
 * This simple library only provides a limited number of operations
 * which are sufficient for our purposes. The underlying matrix library
 * is uBLAS (but the library is designed to be easy to switch to any
 * other underlying library without affecting application code).
 *
 * 1. Each matrix entity is a container of either dense or sparse
 * matrix/vector. 
 *
 * 2. The sparse matrix is a vector of sparse vectors.
 *
 * 3. All matrix are assumed to be row-major.
 * 
 * For all vectors and matrices the underlying data structure is of type
 * unbounded_array
 */

class TheMatrix {
      
public:
      
  // constructor
  TheMatrix();
  TheMatrix(const TheMatrix& tm);
  TheMatrix(const TheMatrix& tm, const int& _dense_or_sparse);
  TheMatrix(const unsigned int & _row, const unsigned int& _column, const int& _dense_or_sparse=SML::DENSE, const unsigned int& non_zeros=0);
  
  // destructor
  virtual ~TheMatrix();
  
  // methods
  int  ID() const;                            // the (unique) id number of this object
  int  Type() const;                          // get the type of internal data structure
  // ostream& write(ostream &os) const;          // print the content of the object
  // istream& read(istream &os);                 // read the contents of the object
  void Print() const;      
  void Zero();                                // Set all elements to zero
  void Shape(unsigned int& row, unsigned int& column) const;    // get the shape of the matrix/vector
  unsigned int  Rows() const;                          // number of rows
  unsigned int  Cols() const;                          // number of columns
  unsigned int  Length() const;                        // get the number of elements in the object  
  double Density() const;                     // ratio of nnz/length: max=1, min=0
  double* Data() const;                       // get the mem location of the elements of the object
  
  void Get(const unsigned int& index, double& value) const;                     // get vector element
  void Get(const unsigned int& row, const unsigned int& column, double& value) const;    // get matrix element
  
  void Set(const unsigned int& index, const double& value);                     // set vector element
  void Set(const unsigned int& row, const unsigned int& column, const double& value);    // set matrix element
  
  // TODO: rename the operation M*v as "Mult" instead of "Dot"
  void Dot(const TheMatrix& other, TheMatrix& output) const;    // dot product of vec/mat or mat/mat
  void Dot(const TheMatrix& other, double& output) const;         // dot product of two vectors!
  
  void TransposeDot(const TheMatrix& other, TheMatrix& output) const;    // dot product of vec^T/mat or mat/mat
  void TransposeDot(const TheMatrix& other, double& output) const;         // dot product of two vectors!
  
	void Norm1(double& value) const;    // L_1 norm for vector and matrix
	double Norm1() const;
	void Norm2(double& value) const;    // L_2 norm for vector and Frobenius (or Hilbert-Schmidt) norm for matrix
	double Norm2() const;
	void NormInf(double& value) const;  // L_\infty norm for vector and matrix
	double NormInf() const;
	void IndexNormInf(unsigned int& index) const; // index of entries with largest absolute value
	unsigned int IndexNormInf() const;
  
  void AddRow(const unsigned int& rowidx, const unsigned int& length, const double* values);                       // Fill a dense row into the matrix (during matrix initialization)
  void AddRow(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);   // Fill a sparse row into the matrix (during matrix initialization)	 
  void SetRow(const unsigned int& rowidx, const unsigned int& length, const double* values);                       // Set a denserow to the matrix (anytime)
  void SetRow(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);   // Set a sparse row to the matrix (anytime)
  void GetRow(const unsigned int& rowidx, unsigned int& length, double* values) const;                             // Get a dense row from the matrix
  void GetRow(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int*indices) const;                // Get a sparse row from the matrix
  
  void Scale(const double& scalar);                             // scale the vector/matrix
  void ScaleAdd(const double& scalar, const TheMatrix& other);  // scale other and add to *this
  
  void Add(const TheMatrix& other);               // Addition of two TheMatrix's
  void Minus(const TheMatrix& other);             // Subtraction of two TheMatrix's
  void ElementWiseMult(const TheMatrix& other);   // Element-wise multiplication of two TheMatrix's
  void Assign(const TheMatrix& other);            // Assign other to this
  
  void SumRows(TheMatrix& result); // Sum along the rows of a matrix
  
  TheMatrix* CreateMatrixRowView(const unsigned int& index);  // return a new TheMatrix "pointing" the index-th row of *this
  void CreateMatrixRowView(const unsigned int& index, TheMatrix* mat_row);  // Assign a reference "pointing" the index-th row of *this to mat_row
  
private:
  bool isReference;            // is *this a reference to another matrix?
  unsigned int row;                     // number of rows in the matrix
  unsigned int column;                  // number of columns in the matrix
  unsigned int length;                  // row*col
  int vector_or_matrix;        // {vector,matrix}
  int dense_or_sparse;         // {dense,sparse}
  int column_or_row;           // {column-major,row-major}
  int type;                    // type of internal data structure {dv,sv,dm,sm}
  int id;                  // the rank of the object
  static unsigned int nInstance;    // keep the number of instances
  
  
  ublas::vector<double>            *dv;     // dense vector
  ublas::compressed_vector<double> *sv;     // sparse vector
  ublas::matrix<double>            *dm;     // dense matrix
  ublas::compressed_matrix<double> *sm;     // sparse matrix
  ublas::matrix_row<ublas::matrix<double> > *dmr;  // dense matrix row (a vector>
  ublas::matrix_row<ublas::compressed_matrix<double> > *smr;  // sparse matrix row (a vector>
  
  // member function pointers
  void (TheMatrix::*fptr_Dot_1)(const TheMatrix& other, double& output) const;
  void (TheMatrix::*fptr_Dot_2)(const TheMatrix& other, TheMatrix& output) const;
  
  void (TheMatrix::*fptr_TransposeDot_1)(const TheMatrix& other, double& output) const;
  void (TheMatrix::*fptr_TransposeDot_2)(const TheMatrix& other, TheMatrix& output) const;
  
  void (TheMatrix::*fptr_Get_1)(const unsigned int& index, double& value) const;
  void (TheMatrix::*fptr_Get_2)(const unsigned int& row, const unsigned int& column, double& value) const;
  
  void (TheMatrix::*fptr_Set_1)(const unsigned int& index, const double& value);
  void (TheMatrix::*fptr_Set_2)(const unsigned int& row, const unsigned int& column, const double& value);
  
  void (TheMatrix::*fptr_AddRow_1)(const unsigned int& rowidx, const unsigned int& length, const double* values);
  void (TheMatrix::*fptr_AddRow_2)(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);
  
  void (TheMatrix::*fptr_SetRow_1)(const unsigned int& rowidx, const unsigned int& length, const double* values);
  void (TheMatrix::*fptr_SetRow_2)(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);
  
  void (TheMatrix::*fptr_GetRow_1)(const unsigned int& rowidx, unsigned int& length, double* values) const;
  void (TheMatrix::*fptr_GetRow_2)(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int* indices) const;
  
  void (TheMatrix::*fptr_ScaleAdd_1)(const double& scalar, const TheMatrix& other);
  void (TheMatrix::*fptr_Add_1)(const TheMatrix& other);
  void (TheMatrix::*fptr_Minus_1)(const TheMatrix& other);
  void (TheMatrix::*fptr_ElementWiseMult_1)(const TheMatrix& other);
  void (TheMatrix::*fptr_Assign_1)(const TheMatrix& other);
  // void (TheMatrix::*fptr_write)(ostream &os) const;
  // void (TheMatrix::*fptr_read)(istream &is);

  
  // methods
  
  // assign pointers
  void assign_pointers();
  void copy_thematrix(const TheMatrix& tm, const int& _dense_or_sparse);
  
  // dot product
  void dv_dot_v(const TheMatrix& other, double& output) const;      // dv dot vector
  void dv_dot_m(const TheMatrix& other, TheMatrix& output) const;   // dv dot matrix
  void sv_dot_v(const TheMatrix& other, double& output) const;      // sv dot vector
  void sv_dot_m(const TheMatrix& other, TheMatrix& output) const;   // sv dot matrix
  void dm_dot(const TheMatrix& other, TheMatrix& output) const;     // dm dot matrix/vector
  void sm_dot(const TheMatrix& other, TheMatrix& output) const;     // sm dot matrix/vector
  void dmr_dot_v(const TheMatrix& other, double& output) const;      // dmr dot vector
  void dmr_dot_m(const TheMatrix& other, TheMatrix& output) const;      // dmr dot matrix
  void smr_dot_v(const TheMatrix& other, double& output) const;      // smr dot vector
  void smr_dot_m(const TheMatrix& other, TheMatrix& output) const;      // smr dot matrix
  
  // transpose then dot product
  void dm_transposedot(const TheMatrix& other, TheMatrix& output) const;     // dm dot matrix/vector
  void sm_transposedot(const TheMatrix& other, TheMatrix& output) const;     // sm dot matrix/vector
  
  // get a value
  void dv_get(const unsigned int& index, double& value) const;
  void dv_get2(const unsigned int& row, const unsigned int& column, double& value) const;
  void sv_get(const unsigned int& index, double& value) const;
  void sv_get2(const unsigned int& row, const unsigned int& column, double& value) const;
  void dm_get(const unsigned int& row, const unsigned int& column, double& value) const;
  void sm_get(const unsigned int& row, const unsigned int& column, double& value) const;
  void dmr_get(const unsigned int& index, double& value) const;
  void smr_get(const unsigned int& index, double& value) const;
  
  // set a value
  void dv_set(const unsigned int& index, const double& value);
  void dv_set2(const unsigned int& row, const unsigned int& column, const double& value);
  void sv_set(const unsigned int& index, const double& value);
  void sv_set2(const unsigned int& row, const unsigned int& column, const double& value);
  void dm_set(const unsigned int& row, const unsigned int& column, const double& value);
  void sm_set(const unsigned int& row, const unsigned int& column, const double& value);
  void dmr_set(const unsigned int& index, const double& value);
  void smr_set(const unsigned int& index, const double& value);
  
  // add a row
  void dm_addrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values);
  void dm_addrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);
  void sm_addrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values);
  void sm_addrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);
  // dense and sparse vector: addrow == setrow  
  void v_addrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values);
  void v_addrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);
  
  // set a row
  // dm_setrow is equivalent to dm_addrow, so, not implemented (although, we can do some raw mem manipulation to make current function faster...)
  void sm_setrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values);
  void sm_setrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices);  
  
  // get a row
  void dm_getrow_dense(const unsigned int& rowidx, unsigned int& length, double* values) const;
  void dm_getrow_sparse(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int* indices) const;
  void sm_getrow_dense(const unsigned int& rowidx, unsigned int& length, double* values) const;
  void sm_getrow_sparse(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int* indices) const; 
  void v_getrow_dense(const unsigned int& rowidx, unsigned int& length, double* values) const;
  void v_getrow_sparse(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int* indices) const;
  
  // scale other and add to *this
  void dv_scaleadd(const double& scalar, const TheMatrix& other);   
  void sv_scaleadd(const double& scalar, const TheMatrix& other);
  void dm_scaleadd(const double& scalar, const TheMatrix& other); 
  void sm_scaleadd(const double& scalar, const TheMatrix& other);
  void dmr_scaleadd(const double& scalar, const TheMatrix& other); 
  void smr_scaleadd(const double& scalar, const TheMatrix& other);
  
  // addition
  void dv_add(const TheMatrix& other);   
  void sv_add(const TheMatrix& other);  
  void dm_add(const TheMatrix& other);  
  void sm_add(const TheMatrix& other);  
  
  // subtraction
  void dv_minus(const TheMatrix& other);   
  void sv_minus(const TheMatrix& other);  
  void dm_minus(const TheMatrix& other);  
  void sm_minus(const TheMatrix& other);  
  
  // element-wise mult
  void dv_elementwisemult(const TheMatrix& other);   
  void sv_elementwisemult(const TheMatrix& other);
  void dm_elementwisemult(const TheMatrix& other);
  void sm_elementwisemult(const TheMatrix& other);
  
  // assignment
  void dv_assign(const TheMatrix& other);   
  void sv_assign(const TheMatrix& other);  
  void dm_assign(const TheMatrix& other);  
  void sm_assign(const TheMatrix& other);    

}; // TheMatrix
   
/*
 *  INLINED FUNCIONS
 */

/**
 * Return the id of underlying matrix. 
 */
inline int TheMatrix::ID() const{
  return id;
}

/**   
 * Return the type of underlying matrix. 
 */
inline int TheMatrix::Type() const{
  return type;
}

/** 
 * Return the shape of underlying matrix. 
 * 
 * @param row [read/write] 
 * @param column [read/write] 
 */
inline void TheMatrix::Shape(unsigned int& row, unsigned int& column) const{
  row = this->row;
  column = this->column;
  return;
}

/**   
 * Return the number of elements of the object
 */
inline unsigned int TheMatrix::Length() const{
  return length;
}

/**   
 * Return number of rows
 */
inline unsigned int TheMatrix::Rows() const{
   return row;
}

/**   
 * Return the number of columns
 */
inline unsigned int TheMatrix::Cols() const{
  return column;
}

/** Return the ratio of nnz / length
 *  Does not work for matrix_row 
 */
inline double TheMatrix::Density() const
{
        using namespace SML;      
                
        // type assertion
        assert(type != DENSE_MATRIX_ROW && type != SPARSE_MATRIX_ROW);
        double nnz = 0.0;
        double *data = 0;
        
        switch(type)
        {
                case DENSE_VECTOR:
                        data = &(dv->data()[0]);
                        for(unsigned int i=0; i<length; i++) if(SML::abs(data[i]) > SML::ZERO_EPS) nnz += 1;
                        break;
                        
                case DENSE_MATRIX:
                        data = &(dm->data()[0]);
                        for(unsigned int i=0; i<length; i++) if(SML::abs(data[i]) > SML::ZERO_EPS) nnz += 1;
                        break;
                                        
                case SPARSE_VECTOR:
                        nnz = sv->nnz();
                        break;
                      
                case SPARSE_MATRIX:
                        nnz = sm->nnz();
                        break;                        
        }
        
        return nnz/length;
}

/** 
 * Return the memory location of the elements of the object. 
 * Does not work for matrix_row 
 */
inline double* TheMatrix::Data() const 
{                
        // svnvish: BUGBUG
        // deprecated function?
        
        // type assertion
        assert(type != SML::DENSE_MATRIX_ROW && type != SML::SPARSE_MATRIX_ROW);
        double* address = 0;
        
        switch(type)
        {
                case SML::DENSE_VECTOR:  address = (&(dv->data()[0])); break;
                case SML::DENSE_MATRIX:  address = (&(dm->data()[0])); break;
                case SML::SPARSE_VECTOR: address = (&(sv->value_data()[0])); break;
                case SML::SPARSE_MATRIX: address = (&(sm->value_data()[0])); break;  
        }
        return address;
}

/**   
 * Zero out all elements. 
 */ 
inline void TheMatrix::Zero(){
  
  switch(type){
  case SML::DENSE_VECTOR:  (*dv).clear(); break;
  case SML::SPARSE_VECTOR: (*sv).clear(); break;
  case SML::DENSE_MATRIX:  (*dm).clear(); break;
  case SML::SPARSE_MATRIX: (*sm).clear(); break;
  }
  return;
}

/** 
 * Get element at given index of a vector.
 * 
 * @param index [read] index to read data from.
 * @param value [write] data value.
 */
inline void TheMatrix::Get(const unsigned int& index, double& value) const{
  (this->*fptr_Get_1)(index, value);
  return;
}

/** 
 * Get element at given row and column of a matrix. 
 * 
 * @param row [read] row to read data from .
 * @param column [read] column to read data from.
 * @param value [write] data value.
 */
inline void TheMatrix::Get(const unsigned int& row, const unsigned int& column, double& value) const{
  (this->*fptr_Get_2)(row, column, value);
  return;
}

/** 
 * Set element at given index of a vector.
 * 
 * @param index [read] index to read data from.
 * @param value [read] data value.
 */
inline void TheMatrix::Set(const unsigned int& index, const double& value){
  (this->*fptr_Set_1)(index, value);
  return;
}

/** 
 * Set element at given row and column of a matrix. 
 * 
 * @param row [read] row to read data from .
 * @param column [read] column to read data from.
 * @param value [read] data value.
 */
inline void TheMatrix::Set(const unsigned int& row, const unsigned int& column, const double& value){
   (this->*fptr_Set_2)(row, column, value);
   return;
}

/** 
 * Return the one-norm of the object. CAUTION: The matrix one-norm is not
 * the same as the one-norm of a vector. Use this function with care. 
 * 
 * @param value [write] one-norm. 
 */
inline void TheMatrix::Norm1(double& value) const{
  
  switch(type) {
  case SML::DENSE_VECTOR:  value = norm_1((*dv)); break;
  case SML::SPARSE_VECTOR: value = norm_1((*sv)); break;
  case SML::DENSE_MATRIX:  value = norm_1((*dm)); break;
  case SML::SPARSE_MATRIX: value = norm_1((*sm)); break;
  case SML::DENSE_MATRIX_ROW:  value = norm_1((*dmr)); break;
  case SML::SPARSE_MATRIX_ROW:  value = norm_1((*smr)); break;
  }
  return;
}

inline double TheMatrix::Norm1() const
{
	double value = 0.0;
	switch(type) {
	case SML::DENSE_VECTOR:  value = norm_1((*dv)); break;
	case SML::SPARSE_VECTOR: value = norm_1((*sv)); break;
	case SML::DENSE_MATRIX:  value = norm_1((*dm)); break;
	case SML::SPARSE_MATRIX: value = norm_1((*sm)); break;
	case SML::DENSE_MATRIX_ROW:  value = norm_1((*dmr)); break;
	case SML::SPARSE_MATRIX_ROW:  value = norm_1((*smr)); break;
	}
	return value;
}


/** 
 * Return the two-norm of the object. In the case of a vector return the
 * vector norm and in the case of a matrix return the Frobenious norm. 
 * 
 * @param value [write] two-norm. 
 */
inline void TheMatrix::Norm2(double& value) const{
  
  switch(type){
  case SML::DENSE_VECTOR:  value = norm_2((*dv)); break;
  case SML::SPARSE_VECTOR: value = norm_2((*sv)); break;
  case SML::DENSE_MATRIX:  value = norm_frobenius((*dm)); break;
  case SML::SPARSE_MATRIX: value = norm_frobenius((*sm)); break;
  case SML::DENSE_MATRIX_ROW:  value = norm_2((*dmr)); break;
  case SML::SPARSE_MATRIX_ROW:  value = norm_2((*smr)); break;
  }
  return;
}

inline double TheMatrix::Norm2() const
{
	double value = 0.0;
	switch(type){
	case SML::DENSE_VECTOR:  value = norm_2((*dv)); break;
	case SML::SPARSE_VECTOR: value = norm_2((*sv)); break;
	case SML::DENSE_MATRIX:  value = norm_frobenius((*dm)); break;
	case SML::SPARSE_MATRIX: value = norm_frobenius((*sm)); break;
	case SML::DENSE_MATRIX_ROW:  value = norm_2((*dmr)); break;
	case SML::SPARSE_MATRIX_ROW:  value = norm_2((*smr)); break;
	}
	return value;
}


/** 
 * Return the infty-norm of the object. CAUTION: The matrix infty-norm is not
 * the same as the infty-norm of a vector. Use this function with care. 
 * 
 * @param value [write] infty-norm. 
 */
 inline void TheMatrix::NormInf(double& value) const{
   
   switch(type){
   case SML::DENSE_VECTOR:  value = norm_inf((*dv)); break;
   case SML::SPARSE_VECTOR: value = norm_inf((*sv)); break;
   case SML::DENSE_MATRIX:  value = norm_inf((*dm)); break;
   case SML::SPARSE_MATRIX: value = norm_inf((*sm)); break;
   case SML::DENSE_MATRIX_ROW:  value = norm_inf((*dmr)); break;
   case SML::SPARSE_MATRIX_ROW:  value = norm_inf((*smr)); break;
   }
   return;
 }


 inline double TheMatrix::NormInf() const
{
	double value = 0.0;
	switch(type){
	case SML::DENSE_VECTOR:  value = norm_inf((*dv)); break;
	case SML::SPARSE_VECTOR: value = norm_inf((*sv)); break;
	case SML::DENSE_MATRIX:  value = norm_inf((*dm)); break;
	case SML::SPARSE_MATRIX: value = norm_inf((*sm)); break;
	case SML::DENSE_MATRIX_ROW:  value = norm_inf((*dmr)); break;
	case SML::SPARSE_MATRIX_ROW:  value = norm_inf((*smr)); break;
	}
	return value;
 }


/** 
 * Return the index-norm of a vector. Assertion error in case of a
 * matrix.  
 * 
 * @param index [write] index-norm. 
 */
inline void TheMatrix::IndexNormInf(unsigned int& index) const{
  
  assert(type != SML::DENSE_MATRIX && type != SML::SPARSE_MATRIX);
  switch(type){
  case SML::DENSE_VECTOR:  index = index_norm_inf((*dv)); break;
  case SML::SPARSE_VECTOR: index = index_norm_inf((*sv)); break;
  case SML::DENSE_MATRIX_ROW:  index = index_norm_inf((*dmr)); break;
  case SML::SPARSE_MATRIX_ROW: index = index_norm_inf((*smr)); break;
  }
  return;
}

inline unsigned int TheMatrix::IndexNormInf() const
{
	unsigned int index = 0;
	assert(type != SML::DENSE_MATRIX && type != SML::SPARSE_MATRIX);
	switch(type){
	case SML::DENSE_VECTOR:  index = index_norm_inf((*dv)); break;
	case SML::SPARSE_VECTOR: index = index_norm_inf((*sv)); break;
	case SML::DENSE_MATRIX_ROW:  index = index_norm_inf((*dmr)); break;
	case SML::SPARSE_MATRIX_ROW: index = index_norm_inf((*smr)); break;
	}
	return index;
}



/** 
 * Add a dense row to a matrix. 
 * 
 * @param rowidx [read] row index.
 * @param length [read] length of the array.
 * @param values [read] array of values. 
 */
inline void TheMatrix::AddRow(const unsigned int& rowidx, 
                              const unsigned int& length, 
                              const double* values){
  (this->*fptr_AddRow_1)(rowidx, length, values);
  return;
}

/** 
 * Add a sparse row to a matrix. 
 * 
 * @param rowidx [read] row index. 
 * @param length [read] length of the array. 
 * @param values [read] array of values. 
 * @param indices [read] non-zero indices. 
 */
inline void TheMatrix::AddRow(const unsigned int& rowidx, 
                              const unsigned int& length, 
                              const double* values, 
                              const unsigned int* indices){
  // svnvish: BUGBUG
  // why do we need length?
  (this->*fptr_AddRow_2)(rowidx, length, values, indices);
  return;
}

/** 
 * Set a dense row of a matrix. 
 * 
 * @param rowidx [read] row index.
 * @param length [read] length of the array.
 * @param values [read] array of values. 
 */
inline void TheMatrix::SetRow(const unsigned int& rowidx, 
                              const unsigned int& length, 
                              const double* values){
   (this->*fptr_SetRow_1)(rowidx,length,values);
   return;
}

/** 
 * Set a sparse row to a matrix. 
 * 
 * @param rowidx [read] row index. 
 * @param length [read] length of the array. 
 * @param values [read] array of values. 
 * @param indices [read] non-zero indices. 
 */
inline void TheMatrix::SetRow(const unsigned int& rowidx, 
                              const unsigned int& length, 
                              const double* values, 
                              const unsigned int* indices){
  //svnvish: BUGBUG
  // why do we need length?
  (this->*fptr_SetRow_2)(rowidx,length,values,indices);
  return;
}

/** 
 * Get a dense row from a matrix. 
 * 
 * @param rowidx [read] row index.
 * @param length [write] length of the array.
 * @param values [write] array of values. 
 */
inline void TheMatrix::GetRow(const unsigned int& rowidx, 
                              unsigned int& length, 
                              double* values) const{
  (this->*fptr_GetRow_1)(rowidx,length,values);
  return;
}


/** 
 * Get a sparse row of a matrix. 
 * 
 * @param rowidx [read] row index. 
 * @param length [write] length of the array. 
 * @param values [write] array of values. 
 * @param indices [write] non-zero indices. 
 */
inline void TheMatrix::GetRow(const unsigned int& rowidx, 
                              unsigned int& length, 
                              double* values, 
                              unsigned int* indices) const{
  (this->*fptr_GetRow_2)(rowidx,length,values,indices);
  return;
}

/** 
 * Scale all values of the matrix/vector by a scalar.
 * 
 * @param sfactor [read] scaling factor. 
 */
inline void TheMatrix::Scale(const double& sfactor)
{

    switch(type)
    {
        case SML::DENSE_VECTOR:  (*dv) *= sfactor; break;
        case SML::SPARSE_VECTOR: (*sv) *= sfactor; break;
        case SML::DENSE_MATRIX:  (*dm) *= sfactor; break;
        case SML::SPARSE_MATRIX: (*sm) *= sfactor; break;
    }
}

/** 
 * Scale all values of the matrix/vector by a scalar and add another
 * matrix.  
 * 
 * @param sfactor [read] scaling factor. 
 * @param other [read] matrix to add. 
 */
inline void TheMatrix::ScaleAdd(const double& sfactor, 
                                const TheMatrix& other){
   (this->*fptr_ScaleAdd_1)(sfactor, other);
   return;
}

/** 
 * Add another matrix. 
 * 
 * @param other [read] matrix to add. 
 */
inline void TheMatrix::Add(const TheMatrix& other){
   (this->*fptr_Add_1)(other);
   return;
}

/** 
 * Subtract another matrix. 
 * 
 * @param other [read] matrix to add. 
 */
inline void TheMatrix::Minus(const TheMatrix& other){
   (this->*fptr_Minus_1)(other);
   return;
}

/** 
 * Hadamard product with another matrix. Dimensions must confirm. 
 * 
 * @param other [read] matrix to elementwise multiply with. 
 */
inline void TheMatrix::ElementWiseMult(const TheMatrix& other){
  (this->*fptr_ElementWiseMult_1)(other);
  return;
}

/** 
 * Assign another matrix to *this. 
 * 
 * @param other [read] matrix to assign. 
 */
inline void TheMatrix::Assign(const TheMatrix& other){
  (this->*fptr_Assign_1)(other);
  return;
}


/*
 *   INLINED PRIVATE FUNCTIONS. 
 */


/**  
 * Get dense vector element using vector index
 */
inline void TheMatrix::dv_get(const unsigned int& index, 
                              double& value) const{
  value = (*dv)[index];
  return;
}

/**
 * Get dense vector element using matrix indices, i.e. either row or column must be 0.
 */
inline void TheMatrix::dv_get2(const unsigned int& row, 
                               const unsigned int& column, 
                               double& value) const{
  assert((column_or_row == SML::ROW && row == 0) or
         (column_or_row == SML::COLUMN && column == 0));
  value = (*dv)[row+column];
  return;
}


/**
 * Get sparse vector element
 */
inline void TheMatrix::sv_get(const unsigned int& index, 
                              double& value) const{
  value = (*sv)[index];
  return;
}

/**
 * Get sparse vector element using matrix indices, i.e. either row or column must be 0.
 */
inline void TheMatrix::sv_get2(const unsigned int& row, 
                               const unsigned int& column, 
                               double& value) const{
  assert((column_or_row == SML::ROW && row == 0) or
         (column_or_row == SML::COLUMN && column == 0));
  value = (*sv)[row+column];
  return;
}


/**
 * Get dense matrix element
 *   [090707]
 */
inline void TheMatrix::dm_get(const unsigned int& row, 
                              const unsigned int& column, 
                              double& value) const{
  // try not to access matrix elements in this way; it's slow!
  value = (*dm)(row,column);
  return;
}


/**
 * Get sparse matrix element
 *   [090707]
 */
inline void TheMatrix::sm_get(const unsigned int& row, 
                              const unsigned int& column, 
                              double& value) const{
  // try not to access matrix elements in this way; it's slow!
  value = (*sm)(row,column);
  return;
}


/**
 * Get dense matrix row element using vector index
 */
inline void TheMatrix::dmr_get(const unsigned int& index, 
                               double& value) const{
  value = (*dmr)[index];
  return;
}

/**
 * Get sparse matrix row element using vector index
 */
inline void TheMatrix::smr_get(const unsigned int& index, 
                               double& value) const{
  value = (*smr)[index];
  return;
}

/**
 * Set dense vector element
 */
inline void TheMatrix::dv_set(const unsigned int& index, 
                              const double& value){
  (*dv)[index] = value;
  return;
}

/**
 * Set dense vector element using matrix indices, i.e. either row or column must be 0.
 */
inline void TheMatrix::dv_set2(const unsigned int& row, 
                               const unsigned int& column, 
                               const double& value){
  assert((column_or_row == SML::ROW && row == 0) or
         (column_or_row == SML::COLUMN && column == 0));
  (*dv)[row+column] = value;
  return;
}

/**
 * Set sparse vector element
 */
inline void TheMatrix::sv_set(const unsigned int& index, 
                              const double& value){
  (*sv)[index] = value;
  return;
}

/**
 * Set sparse vector element using matrix indices, i.e. either row or column must be 0.
 */
inline void TheMatrix::sv_set2(const unsigned int& row, 
                               const unsigned int& column, 
                               const double& value){
  assert((column_or_row == SML::ROW && row == 0) or
         (column_or_row == SML::COLUMN && column == 0));
  (*sv)[row+column] = value;
  return;
}

/**
 * Set dense matrix element
 *   [090707]
 */
inline void TheMatrix::dm_set(const unsigned int& row, 
                              const unsigned int& column, 
                              const double& value){
  // try not to access matrix elements in this way; it's slow!
  (*dm)(row,column) = value;
  return;
}

/**
 * Set sparse matrix element
 *   [090707]
 */
inline void TheMatrix::sm_set(const unsigned int& row, 
                              const unsigned int& column, 
                              const double& value){
  // try not to access matrix elements in this way; it's slow!
  (*sm)(row,column) = value;
  return;
}

/**
 * Set dense matrix row element
 */
inline void TheMatrix::dmr_set(const unsigned int& index, 
                               const double& value){
  (*dmr)[index] = value;
  return;
}

/**
 * Set sparse matrix row element
 */
inline void TheMatrix::smr_set(const unsigned int& index, 
                               const double& value){
  (*smr)[index] = value;
  return;
}

/**
 * Scale other and add to dense vector
 */
inline void TheMatrix::dv_scaleadd(const double& scalar, 
                                   const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::DENSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  (*dv) += scalar * (*other.dv); break;
  case SML::SPARSE_VECTOR: (*dv) += scalar * (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: (*dv) += scalar * (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: (*dv) += scalar * (*other.smr); break;
  }
  return;
}

/**
 * Scale other and add to dense matrix row view
 */
inline void TheMatrix::dmr_scaleadd(const double& scalar, const TheMatrix& other)
{
        int otherType = other.Type();
        assert(id != other.ID());      // use Scale() instead!
        assert(type == SML::DENSE_MATRIX_ROW);
        assert(otherType == SML::DENSE_VECTOR ||
              otherType == SML::SPARSE_VECTOR ||
              otherType == SML::DENSE_MATRIX_ROW ||
              otherType == SML::SPARSE_MATRIX_ROW);
        switch(otherType)
        {
                case SML::DENSE_VECTOR:  noalias(*dmr) += scalar * (*other.dv); break;
                case SML::SPARSE_VECTOR: noalias(*dmr) += scalar * (*other.sv); break;
                case SML::DENSE_MATRIX_ROW: noalias(*dmr) += scalar * (*other.dmr); break;
                case SML::SPARSE_MATRIX_ROW: noalias(*dmr) += scalar * (*other.smr); break;
        }
}


/**
 * Add other to dense vector
 */
inline void TheMatrix::dv_add(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::DENSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  (*dv) += (*other.dv); break;
  case SML::SPARSE_VECTOR: (*dv) += (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: (*dv) += (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: (*dv) += (*other.smr); break;
  }
  return;
}

/**
 * Subtract other from dense vector
 */
inline void TheMatrix::dv_minus(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      
  assert(type == SML::DENSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  (*dv) -= (*other.dv); break;
  case SML::SPARSE_VECTOR: (*dv) -= (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: (*dv) -= (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: (*dv) -= (*other.smr); break;
  }
  return;
}

/**
 * Element-wise multiply other with dense vector
 */
inline void TheMatrix::dv_elementwisemult(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      
  assert(type == SML::DENSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  (*dv) = element_prod(*dv, *(other.dv)); break;
  case SML::SPARSE_VECTOR: (*dv) = element_prod(*dv, *(other.sv)); break;
  case SML::DENSE_MATRIX_ROW: (*dv) = element_prod(*dv, *(other.dmr)); break;
  case SML::SPARSE_MATRIX_ROW: (*dv) = element_prod(*dv, *(other.smr)); break;
  }
  return;
}

/**
 * Scale other and add to sparse vector
 */
inline void TheMatrix::sv_scaleadd(const double& scalar, 
                                   const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::SPARSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  noalias(*sv) += scalar * (*other.dv); break;
  case SML::SPARSE_VECTOR: noalias(*sv) += scalar * (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: noalias(*sv) += scalar * (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: noalias(*sv) += scalar * (*other.smr); break;
  }
  return;
}


/**
 * Scale other and add to sparse matrix row view
 */
inline void TheMatrix::smr_scaleadd(const double& scalar, 
                                   const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::SPARSE_MATRIX_ROW);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  noalias(*smr) += scalar * (*other.dv); break;
  case SML::SPARSE_VECTOR: noalias(*smr) += scalar * (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: noalias(*smr) += scalar * (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: noalias(*smr) += scalar * (*other.smr); break;
  }
  return;
}


/**
 * Add other to sparse vector
 */
inline void TheMatrix::sv_add(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::SPARSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  (*sv) += (*other.dv); break;
  case SML::SPARSE_VECTOR: (*sv) += (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: (*sv) += (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: (*sv) += (*other.smr); break;
  }
  return;
}

/**
 * Subtract other from sparse vector
 */
inline void TheMatrix::sv_minus(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());
  assert(type == SML::SPARSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  
  switch(otherType){
  case SML::DENSE_VECTOR:  (*sv) -= (*other.dv); break;
  case SML::SPARSE_VECTOR: (*sv) -= (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: (*sv) -= (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: (*sv) -= (*other.smr); break;
  }
  return;
}

/**
 * Element-wise multiply other with sparse vector
 */
inline void TheMatrix::sv_elementwisemult(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      
  assert(type == SML::SPARSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR or
         otherType == SML::SPARSE_VECTOR or
         otherType == SML::DENSE_MATRIX_ROW or
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  (*sv) = element_prod(*sv, *other.dv); break;
  case SML::SPARSE_VECTOR: (*sv) = element_prod(*sv, *other.sv); break;
  case SML::DENSE_MATRIX_ROW: (*sv) = element_prod(*sv, *other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: (*sv) = element_prod(*sv, *other.smr); break;
  }
  return;
}


/**
 * Scale other and add to dense matrix
 */
inline void TheMatrix::dm_scaleadd(const double& scalar, 
                                   const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::DENSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX || 
         otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  noalias((*dm)) += scalar * (*other.dm); break;
  case SML::SPARSE_MATRIX: noalias((*dm)) += scalar * (*other.sm); break;
  }
  return;
}

/**
 * Add other to dense matrix
 */
inline void TheMatrix::dm_add(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::DENSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX || 
         otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  (*dm) += (*other.dm); break;
  case SML::SPARSE_MATRIX: (*dm) += (*other.sm); break;
  }
  return;
}

/**
 * Subtract other from dense matrix
 */
inline void TheMatrix::dm_minus(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());
  assert(type == SML::DENSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX || 
         otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  (*dm) -= (*other.dm); break;
  case SML::SPARSE_MATRIX: (*dm) -= (*other.sm); break;
  }
  return;
}

/**
 * Element-wise multiply other with dense matrix
 */
inline void TheMatrix::dm_elementwisemult(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      
  assert(type == SML::DENSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX || 
         otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  (*dm) = element_prod(*dm, *other.dm); break;
  case SML::SPARSE_MATRIX: (*dm) = element_prod(*dm, *other.sm); break;
  }
  return;
}


/**
 * Scale other and add to sparse matrix
 */
inline void TheMatrix::sm_scaleadd(const double& scalar, 
                                   const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::SPARSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX || 
         otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  noalias((*sm)) += scalar * (*other.dm); break;
  case SML::SPARSE_MATRIX: noalias((*sm)) += scalar * (*other.sm); break;
  }
  return;
}

/**
 * Add other to sparse matrix
 */
inline void TheMatrix::sm_add(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      // use Scale() instead!
  assert(type == SML::SPARSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX  || 
         otherType == SML::SPARSE_MATRIX);
  
  switch(otherType){
  case SML::DENSE_MATRIX:  (*sm) += (*other.dm); break;
  case SML::SPARSE_MATRIX: (*sm) += (*other.sm); break;
  }
  return;
}

/**
 * Subtract other from sparse matrix
 */
inline void TheMatrix::sm_minus(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());
  assert(type == SML::SPARSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX || 
         otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  (*sm) -= (*other.dm); break;
  case SML::SPARSE_MATRIX: (*sm) -= (*other.sm); break;
  }
  return;
}

/**
 * Element-wise multiply other with sparse matrix
 */
inline void TheMatrix::sm_elementwisemult(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      
  assert(type == SML::SPARSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX  || 
         otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  (*sm) = element_prod(*sm, *other.dm); break;
  case SML::SPARSE_MATRIX: (*sm) = element_prod(*sm, *other.sm); break;
  }
  return;
}


/**
 * Assign other to dense vector
 *    other could be a dense matrix, dense vector or sparse vector.
 *    For the case of dense matrix, memcpy is used
 */
inline void TheMatrix::dv_assign(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());     
  assert(type == SML::DENSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR  ||
         otherType == SML::SPARSE_VECTOR || 
         otherType == SML::DENSE_MATRIX_ROW || 
         otherType == SML::SPARSE_MATRIX_ROW or
         otherType == SML::DENSE_MATRIX);
  
  switch(otherType){
  case SML::DENSE_VECTOR:  noalias(*dv) = (*other.dv); break;
  case SML::SPARSE_VECTOR: noalias(*dv) = (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: noalias(*dv) = (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: noalias(*dv) = (*other.smr); break;
    
  case SML::DENSE_MATRIX: 
    assert(length == other.length);
    memcpy(&(dv->data()[0]), &(other.dm->data()[0]), sizeof(double)*length);
    break;
  }
  return;
}


/**
 * Assign other to sparse vector
 */
inline void TheMatrix::sv_assign(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());     
  assert(type == SML::SPARSE_VECTOR);
  assert(otherType == SML::DENSE_VECTOR  || 
         otherType == SML::SPARSE_VECTOR || 
         otherType == SML::DENSE_MATRIX_ROW || 
         otherType == SML::SPARSE_MATRIX_ROW);
  switch(otherType){
  case SML::DENSE_VECTOR:  noalias(*sv) = (*other.dv); break;
  case SML::SPARSE_VECTOR: noalias(*sv) = (*other.sv); break;
  case SML::DENSE_MATRIX_ROW: noalias(*sv) = (*other.dmr); break;
  case SML::SPARSE_MATRIX_ROW: noalias(*sv) = (*other.smr); break;
  }
  return;
}


/**
 * Assign other to dense matrix
 * other can be a dense vector, dense matrix, or sparse matrix
 * for the case of dense vector, memcpy is used.
 */
inline void TheMatrix::dm_assign(const TheMatrix& other){
   
   int otherType = other.Type();
   assert(id != other.ID());      
   assert(type == SML::DENSE_MATRIX);
   assert(otherType == SML::DENSE_MATRIX  or
          otherType == SML::SPARSE_MATRIX or
          otherType == SML::DENSE_VECTOR);
   switch(otherType){
      case SML::DENSE_MATRIX:  noalias(*dm) = (*other.dm); break;
      case SML::SPARSE_MATRIX: noalias(*dm) = (*other.sm); break;
      case SML::DENSE_VECTOR: 
         assert(length == other.length);
         memcpy(&(dm->data()[0]), &(other.dv->data()[0]), sizeof(double)*length);
         break;
   }
   return;
}


/**   
 * Assign other to sparse matrix
 */
inline void TheMatrix::sm_assign(const TheMatrix& other){
  int otherType = other.Type();
  assert(id != other.ID());      
  assert(type == SML::SPARSE_MATRIX);
  assert(otherType == SML::DENSE_MATRIX  || otherType == SML::SPARSE_MATRIX);
  switch(otherType){
  case SML::DENSE_MATRIX:  noalias(*sm) = (*other.dm); break;
  case SML::SPARSE_MATRIX: noalias(*sm) = (*other.sm); break;
  }
  return;
}

void printSparseVector(unsigned int dim, TheMatrix * v);

/*
  How to add a new method:
  1. Declare a public method name
  2. Declare and define private functions tailored to different combination of types of TheMatrix
  3. Declare member function pointer of type (2)
  4. Direct the member function pointer to "correct" private function defined in (2) in constructor
  5. Call the member function pointer in the public method declared in (1)
*/
#endif

