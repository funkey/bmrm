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
 * Last Updated: (15/01/2008)   
 */


/*   Simple Matrix Library (SML)
 *
 *   1. Each matrix entity is a container of either dense or sparse matrix/vector.
 *   2. The sparse matrix is a vector of sparse vectors.
 *   3. All matrices are assumed to be row-major,  
 */

#ifndef _SML_CPP_
#define _SML_CPP_

#include <iostream>
#include "sml.hpp"


using namespace std;
using namespace SML;
namespace ublas = boost::numeric::ublas;



// who added this?
void printSparseVector(unsigned int dim, TheMatrix * v)
{
	for(unsigned int i=0;i<dim;i++)
	{
		double val;
		v->Get(0,i,val);
		if(fabs(val) > 1e-20)
		{
			cout<<i<<" : "<<val<<"\t";
		}
	}
}


/**  Constructor
 */
TheMatrix::TheMatrix()
   : isReference(false), 
     type(0),
     dv(0),
     sv(0),
     dm(0),
     sm(0),
     dmr(0),
     smr(0)
{
        id = nInstance;
        nInstance++;		
#ifdef DEBUG
        std::cout << "in constructor 1:  TheMatrix ID: " << id << std::endl;
#endif

}


/**  Copy constructor
 */
TheMatrix::TheMatrix(const TheMatrix& tm)
   : isReference(false), 
     type(0),
     dv(0),
     sv(0),
     dm(0),
     sm(0),
     dmr(0),
     smr(0)
{
        copy_thematrix(tm, tm.dense_or_sparse);
       
        id = nInstance;
        nInstance++;
        
#ifdef DEBUG
        std::cout << "in constructor 2:  TheMatrix ID: " << id << std::endl;
        cout.flush();
#endif

}

/**  Copy constructor
 *   with option to convert the type (e.g. dense or sparse) of the object into another
 */
TheMatrix::TheMatrix(const TheMatrix& tm, const int& _dense_or_sparse)
   : isReference(false),
     type(0),
     dv(0),
     sv(0),
     dm(0),
     sm(0),
     dmr(0),
     smr(0)
{
        copy_thematrix(tm, _dense_or_sparse);

        id = nInstance;
        nInstance++;

#ifdef DEBUG
        std::cout << "in constructor 3: TheMatrix ID: " << id << std::endl;
        cout.flush();
#endif

}


/**   Copy create a new thematrix with the same content as tm except the type.
 */
void TheMatrix::copy_thematrix(const TheMatrix& tm, const int& _dense_or_sparse)
{
        row = tm.row;
        column = tm.column;
        length = tm.length;
   
        vector_or_matrix = tm.vector_or_matrix;
        dense_or_sparse = _dense_or_sparse;
   
        assert(row > 0 && column > 0);
        if(column == 1 || row == 1) 
        {
                vector_or_matrix = SML::VECTOR;
                if(row == 1) 
                        column_or_row = SML::ROW;
                else
                        column_or_row = SML::COLUMN;
        }
   
        // determine the type
        switch(vector_or_matrix) 
        {
                case SML::VECTOR: 
                        switch(dense_or_sparse) 
                        {
                                case SML::DENSE: type = SML::DENSE_VECTOR; break;
                                case SML::SPARSE:  type = SML::SPARSE_VECTOR; break;
                        }
                        break;
         
                case SML::MATRIX:
                        switch(dense_or_sparse) 
                        {
                                case SML::DENSE: type = SML::DENSE_MATRIX; break;
                                case SML::SPARSE: type = SML::SPARSE_MATRIX; break;
                        }
                        break;
        }
   
        // allocate new ublas object with initialization values
        switch(tm.type) 
        {
                case SML::DENSE_VECTOR: 
                        switch(type) 
                        {
                                case SML::DENSE_VECTOR: dv = new ublas::vector<double>(*(tm.dv)); break;
                                case SML::SPARSE_VECTOR: sv = new ublas::compressed_vector<double>(*(tm.dv)); break;
                        }
                        break;
         
                case SML::SPARSE_VECTOR:
                        switch(type) 
                        {
                                case SML::DENSE_VECTOR: dv = new ublas::vector<double>(*(tm.sv)); break;
                                case SML::SPARSE_VECTOR: sv = new ublas::compressed_vector<double>(*(tm.sv)); break;
                        }
                        break;

                case SML::DENSE_MATRIX: 
                        switch(type) 
                        {
                                case SML::DENSE_MATRIX: dm = new ublas::matrix<double>(*(tm.dm)); break;
                                case SML::SPARSE_MATRIX: sm = new ublas::compressed_matrix<double>(*(tm.dm)); break;
                        }
                        break;
      
                case SML::SPARSE_MATRIX: 
                        switch(type) 
                        {
                                case SML::DENSE_MATRIX: dm = new ublas::matrix<double>(*(tm.sm)); break;
                                case SML::SPARSE_MATRIX: sm = new ublas::compressed_matrix<double>(*(tm.sm)); break;
                        }
                        break;

                case SML::DENSE_MATRIX_ROW: 
                        switch(type) 
                        {
                                case SML::DENSE_VECTOR: dv = new ublas::vector<double>(*(tm.dmr)); break;
                                case SML::SPARSE_VECTOR: sv = new ublas::compressed_vector<double>(*(tm.dmr)); break;
                        }
                        break;
      
                case SML::SPARSE_MATRIX_ROW: 
                        switch(type) 
                        {
                                case SML::DENSE_VECTOR: dv = new ublas::vector<double>(*(tm.smr)); break;
                                case SML::SPARSE_VECTOR: sv = new ublas::compressed_vector<double>(*(tm.smr)); break;
                        }
                        break;

        }
   
        assign_pointers();
}


/** Constructor.
 *  1. If either row or col is 1, the object will be a vector.
 *  2. All vector/matrix elements will be initialized to zeros.
 */
TheMatrix::TheMatrix(const unsigned int & _row, const unsigned int& _column, const int& _dense_or_sparse, const unsigned int& non_zeros)
   : isReference(false),
     row(_row),
     column(_column),
     length(row*column),
     vector_or_matrix(SML::MATRIX),
     dense_or_sparse(_dense_or_sparse),
     column_or_row(ROW),
     type(0),
     dv(0),
     sv(0),
     dm(0),
     sm(0),
     dmr(0),
     smr(0)
{
        assert(row > 0 && column > 0);
        if(column == 1 || row == 1) 
        {
                vector_or_matrix = SML::VECTOR;
                if(row == 1) 
                        column_or_row = SML::ROW;
                else
                        column_or_row = SML::COLUMN;
        }

        switch(vector_or_matrix)
        {
                case SML::VECTOR: 
                        switch(dense_or_sparse) 
                        {
                                case SML::DENSE:
                                        dv = new ublas::vector<double>(length);
                                        dv->clear();
                                        type = SML::DENSE_VECTOR;  
                                        break;
               
                                case SML::SPARSE:
                                        sv = new ublas::compressed_vector<double>(length,non_zeros);
                                        sv->clear();
                                        type = SML::SPARSE_VECTOR;  
                                        break;
                        }
                        break;
         
                case SML::MATRIX:
                        switch(dense_or_sparse) 
                        {
                                case SML::DENSE:
                                        dm = new ublas::matrix<double>(row,column);
                                        dm->clear();
                                        type = SML::DENSE_MATRIX;
                                        break;
               
                                case SML::SPARSE:
                                        sm = new ublas::compressed_matrix<double>(row,column,non_zeros);
                                        sm->clear();
                                        type = SML::SPARSE_MATRIX;
                                        break;
                        }
                        break;
        }
        
        assign_pointers();
        id = nInstance;
        nInstance++;   
        
#ifdef DEBUG
        std::cout << "in constructor 4. TheMatrix ID: " << id << std::endl;
        cout.flush();
#endif
} // TheMatrix::TheMatrix


/**   Assign function pointers properly. Called by constructors only.
 */
void TheMatrix::assign_pointers()
{
        switch(type)
        {
                case SML::DENSE_VECTOR: 
                        fptr_Get_1 = &TheMatrix::dv_get;
                        fptr_Get_2 = &TheMatrix::dv_get2;
                        fptr_Set_1 = &TheMatrix::dv_set;   
                        fptr_Set_2 = &TheMatrix::dv_set2;       
                        fptr_AddRow_1 = &TheMatrix::v_addrow_dense;
                        fptr_AddRow_2 = &TheMatrix::v_addrow_sparse;
                        fptr_SetRow_1 = &TheMatrix::v_addrow_dense;
                        fptr_SetRow_2 = &TheMatrix::v_addrow_sparse;
                        fptr_GetRow_1 = &TheMatrix::v_getrow_dense;
                        fptr_GetRow_2 = &TheMatrix::v_getrow_sparse;
                        fptr_TransposeDot_1 = &TheMatrix::dv_dot_v;
                        fptr_Dot_1 = &TheMatrix::dv_dot_v;
                        fptr_Dot_2 = &TheMatrix::dv_dot_m;
                        fptr_TransposeDot_2 = &TheMatrix::dv_dot_m;
                        fptr_ScaleAdd_1 = &TheMatrix::dv_scaleadd;
                        fptr_Add_1 = &TheMatrix::dv_add;
                        fptr_Minus_1 = &TheMatrix::dv_minus;
                        fptr_ElementWiseMult_1 = &TheMatrix::dv_elementwisemult;
                        fptr_Assign_1 = &TheMatrix::dv_assign;
                        // fptr_write = &TheMatrix::dv_write;
                        // fptr_read = &TheMatrix::dv_read;
                        break;
    
                case SML::SPARSE_VECTOR:
                        fptr_Get_1 = &TheMatrix::sv_get;  
                        fptr_Get_2 = &TheMatrix::sv_get2;  
                        fptr_Set_1 = &TheMatrix::sv_set; 
                        fptr_Set_2 = &TheMatrix::sv_set2;
                        fptr_AddRow_1 = &TheMatrix::v_addrow_dense;
                        fptr_AddRow_2 = &TheMatrix::v_addrow_sparse;
                        fptr_SetRow_1 = &TheMatrix::v_addrow_dense;
                        fptr_SetRow_2 = &TheMatrix::v_addrow_sparse;
                        fptr_GetRow_1 = &TheMatrix::v_getrow_dense;
                        fptr_GetRow_2 = &TheMatrix::v_getrow_sparse;
                        fptr_Dot_1 = &TheMatrix::sv_dot_v;
                        fptr_TransposeDot_1 = &TheMatrix::sv_dot_v;
                        fptr_Dot_2 = &TheMatrix::sv_dot_m;
                        fptr_TransposeDot_2 = &TheMatrix::sv_dot_m;
                        fptr_ScaleAdd_1 = &TheMatrix::sv_scaleadd;
                        fptr_Add_1 = &TheMatrix::sv_add;
                        fptr_Minus_1 = &TheMatrix::sv_minus;
                        fptr_ElementWiseMult_1 = &TheMatrix::sv_elementwisemult;
                        fptr_Assign_1 = &TheMatrix::sv_assign;
                        // fptr_write = &TheMatrix::sv_write;
                        // fptr_read = &TheMatrix::sv_read;
                        break;
    
                case SML::DENSE_MATRIX:
                        fptr_Get_2    = &TheMatrix::dm_get;
                        fptr_Set_2    = &TheMatrix::dm_set;
                        fptr_Dot_2    = &TheMatrix::dm_dot;
                        fptr_TransposeDot_2 = &TheMatrix::dm_transposedot;
                        fptr_AddRow_1 = &TheMatrix::dm_addrow_dense;
                        fptr_AddRow_2 = &TheMatrix::dm_addrow_sparse;
                        fptr_SetRow_1 = &TheMatrix::dm_addrow_dense;   // dm_addrow and dm_setrow happen to be identical ;)
                        fptr_SetRow_2 = &TheMatrix::dm_addrow_sparse;  // dm_addrow and dm_setrow happen to be identical ;)
                        fptr_GetRow_1 = &TheMatrix::dm_getrow_dense;
                        fptr_GetRow_2 = &TheMatrix::dm_getrow_sparse;
                        fptr_ScaleAdd_1 = &TheMatrix::dm_scaleadd;
                        fptr_Add_1 = &TheMatrix::dm_add;
                        fptr_Minus_1 = &TheMatrix::dm_minus;
                        fptr_ElementWiseMult_1 = &TheMatrix::dm_elementwisemult;
                        fptr_Assign_1 = &TheMatrix::dm_assign;
                        // fptr_write = &TheMatrix::dm_write;
                        // fptr_read = &TheMatrix::dm_read;
                        break;
    
                case SML::SPARSE_MATRIX:
                        fptr_Get_2    = &TheMatrix::sm_get;
                        fptr_Set_2    = &TheMatrix::sm_set;
                        fptr_Dot_2    = &TheMatrix::sm_dot;
                        fptr_TransposeDot_2 = &TheMatrix::sm_transposedot;
                        fptr_AddRow_1 = &TheMatrix::sm_addrow_dense;
                        fptr_AddRow_2 = &TheMatrix::sm_addrow_sparse;
                        fptr_SetRow_1 = &TheMatrix::sm_setrow_dense;
                        fptr_SetRow_2 = &TheMatrix::sm_setrow_sparse;
                        fptr_GetRow_1 = &TheMatrix::sm_getrow_dense;
                        fptr_GetRow_2 = &TheMatrix::sm_getrow_sparse;
                        fptr_ScaleAdd_1 = &TheMatrix::sm_scaleadd;
                        fptr_Add_1 = &TheMatrix::sm_add;
                        fptr_Minus_1 = &TheMatrix::sm_minus;
                        fptr_ElementWiseMult_1 = &TheMatrix::sm_elementwisemult;
                        fptr_Assign_1 = &TheMatrix::sm_assign;
                        // fptr_write = &TheMatrix::sm_write;
                        // fptr_read = &TheMatrix::sm_read;
                        break;
    
                case SML::DENSE_MATRIX_ROW:
                        fptr_Get_1 = &TheMatrix::dmr_get;
                        fptr_Get_2 = 0;
                        fptr_Set_1 = &TheMatrix::dmr_set;   
                        fptr_Set_2 = 0;
                        fptr_Dot_1 = &TheMatrix::dmr_dot_v;
                        fptr_TransposeDot_1 = &TheMatrix::dmr_dot_v;
                        fptr_Dot_2 = &TheMatrix::dmr_dot_m;
                        fptr_TransposeDot_2 = &TheMatrix::dmr_dot_m;
                        fptr_ScaleAdd_1 = &TheMatrix::dmr_scaleadd;
                        fptr_Add_1 = 0;
                        fptr_Minus_1 = 0;
                        fptr_ElementWiseMult_1 = 0;
                        fptr_Assign_1 = 0;
                        // fptr_write = &TheMatrix::dmr_write;
                        // fptr_read = &TheMatrix::dmr_read;
                        break;
    
                case SML::SPARSE_MATRIX_ROW:
                        fptr_Get_1 = &TheMatrix::smr_get;
                        fptr_Get_2 = 0;
                        fptr_Set_1 = &TheMatrix::smr_set;   
                        fptr_Set_2 = 0;
                        fptr_Dot_1 = &TheMatrix::smr_dot_v;
                        fptr_TransposeDot_1 = &TheMatrix::smr_dot_v;
                        fptr_Dot_2 = &TheMatrix::smr_dot_m;
                        fptr_TransposeDot_2 = &TheMatrix::smr_dot_m;
                        fptr_ScaleAdd_1 = &TheMatrix::smr_scaleadd;
                        fptr_Add_1 = 0;
                        fptr_Minus_1 = 0;
                        fptr_ElementWiseMult_1 = 0;
                        fptr_Assign_1 = 0;
                        // fptr_write = &TheMatrix::smr_write;
                        // fptr_read = &TheMatrix::smr_read;
                        break;
        }     
}


/** Destructor
 */
TheMatrix::~TheMatrix()
{
        if(dv) delete dv;
        if(sv) delete sv;
        if(dm) delete dm;
        if(sm) delete sm;
        if(dmr) delete dmr;
        if(smr) delete smr;
  
        //nInstance--;    // one instance destructed

#ifdef DEBUG
        std::cout << "In destructor: TheMatrix ID: " << id << std::endl;
        cout.flush();
#endif        
  
} // TheMatrix::~TheMatrix



/*
 *   STATIC CLASS VARIABLE
 */
unsigned int TheMatrix::nInstance = 0;



/**   A visible shell to actual dot product (with a vector).
 *    This function should only be called by a vector object.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::Dot(const TheMatrix& other, double& output) const
{
        // pointer to member function
        (this->*fptr_Dot_1)(other, output);
}

/**   A visible shell to actual dot product (with a matrix).
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::Dot(const TheMatrix& other, TheMatrix& output) const
{
        // pointer to member function
        (this->*fptr_Dot_2)(other, output);
}


/**   A visible shell to actual transpose-then-dot product (with a vector).
 *    This function should only be called by a vector object.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::TransposeDot(const TheMatrix& other, double& output) const
{
        // pointer to member function
        (this->*fptr_TransposeDot_1)(other, output);
}

/**   A visible shell to actual transpose-then-dot product (with a matrix).
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::TransposeDot(const TheMatrix& other, TheMatrix& output) const
{
        // pointer to member function
        (this->*fptr_TransposeDot_2)(other, output);
}


/**   Print the content of the object
 *  
 */
void TheMatrix::Print() const
{
        switch(type)
        {
                case SML::DENSE_VECTOR: std::cout << (*dv) << std::endl; break;
                case SML::DENSE_MATRIX: std::cout << (*dm) << std::endl; break;
                case SML::SPARSE_VECTOR: std::cout << (*sv) << std::endl; break;
                case SML::SPARSE_MATRIX: std::cout << (*sm) << std::endl; break;
                case SML::DENSE_MATRIX_ROW: std::cout << (*dmr) << std::endl; break;
                case SML::SPARSE_MATRIX_ROW: std::cout << (*smr) << std::endl; break;
        }
}

/** Sum along the rows of a matrix.(A syntactic sugar)
 *
 *  \param output [write] The output vector
 */ 
void TheMatrix::SumRows(TheMatrix& output)
{
        assert(output.Type() == SML::DENSE_VECTOR); 
        assert(output.Length() == column);
        
        switch(type)
        {
                case SML::DENSE_VECTOR:
                case SML::SPARSE_VECTOR:
                case SML::DENSE_MATRIX_ROW:
                case SML::SPARSE_MATRIX_ROW:
                        output.Assign(*this);
        }
          
        // svnvish: BUGBUG
        // Maybe boost already has such a utility function?
        TheMatrix one(column, 1, SML::DENSE);
        for(unsigned int i = 0; i < column; i++)
        {
                one.Set(i, 1);
        }
        Dot(one, output);
  
        //chteo: [301207:2247] use identity_vector<double> one instead of TheMatrix
}


/*
 *      "ACTUAL" FUNCTIONS HIDDEN FROM USERS
 */


/**   Actual dot product between dv and any dense/sparse vector.
 *    The result is a double value.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::dv_dot_v(const TheMatrix& other, double& output) const
{
    int otherType = other.Type();
    assert(otherType == SML::DENSE_VECTOR ||
           otherType == SML::SPARSE_VECTOR ||
           otherType == SML::DENSE_MATRIX_ROW ||
           otherType == SML::SPARSE_MATRIX_ROW);
    
    switch(otherType) 
    {
    case SML::DENSE_VECTOR:  output = inner_prod(*dv, *(other.dv)); break;
    case SML::SPARSE_VECTOR: output = inner_prod(*dv, *(other.sv)); break;
    case SML::DENSE_MATRIX_ROW: output = inner_prod(*dv, *(other.dmr)); break;
    case SML::SPARSE_MATRIX_ROW: output = inner_prod(*dv, *(other.smr)); break;
    }    
}


/**   Actual dot product between dv and any dense/sparse matrix.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::dv_dot_m(const TheMatrix& other, TheMatrix& output) const
{
        int otherType  = other.Type();
        int outputType = output.Type();
   
        // type assertion
        assert(otherType == SML::DENSE_MATRIX ||  
               otherType == SML::SPARSE_MATRIX || 
               otherType == SML::DENSE_VECTOR ||
               otherType == SML::SPARSE_VECTOR ||
               otherType == SML::DENSE_MATRIX_ROW ||
               otherType == SML::SPARSE_MATRIX_ROW);
        
        assert(outputType == SML::DENSE_VECTOR ||
               outputType == SML::SPARSE_VECTOR ||
               outputType == SML::DENSE_MATRIX_ROW ||
               outputType == SML::SPARSE_MATRIX_ROW);
   
        switch(otherType)
        {
                case SML::SPARSE_MATRIX: 
                        switch(outputType)
                        {
                                case SML::DENSE_VECTOR: axpy_prod(*dv, *(other.sm), *(output.dv)); break;
                                //axpy_prod(ublas::trans(*other.sm), *dv, *(output.dv)); break; //chteo: this is slower!
                                case SML::SPARSE_VECTOR: axpy_prod(trans(*(other.sm)), *dv, *(output.sv)); break;
                                //case SML::SPARSE_VECTOR: noalias(*(output.sv)) = prod(*dv, *(other.sm)); break;
                        }  
                        break;
         
                case SML::DENSE_MATRIX:
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*dv, *(other.dm), *(output.dv)); break;      
                                case SML::SPARSE_VECTOR: axpy_prod(*dv, *(other.dm), *(output.sv)); break;      
                        }
                        break;
         
                // the following cases are for pathological case: output is a 1x1 matrix
                case SML::DENSE_VECTOR:
                case SML::SPARSE_VECTOR:      
                case SML::DENSE_MATRIX_ROW:
                case SML::SPARSE_MATRIX_ROW:             
                        double value = 0.0;
                        dv_dot_v(other, value);
                        switch(outputType)
                        {
                                case SML::DENSE_VECTOR: (*output.dv)[0] = value; break;
                                case SML::SPARSE_VECTOR: (*output.sv)[0] = value; break;
                                case SML::DENSE_MATRIX_ROW: (*output.dmr)[0] = value; break;
                                case SML::SPARSE_MATRIX_ROW: (*output.smr)[0] = value; break;                      
                        }      
                        break;
        }
}

/**   Actual dot product between sv and any dense/sparse vector.
 *    The result is a double value.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::sv_dot_v(const TheMatrix& other, double& output) const
{
    int otherType = other.Type();
   
    assert(otherType == SML::DENSE_VECTOR ||
           otherType == SML::SPARSE_VECTOR ||
           otherType == SML::DENSE_MATRIX_ROW ||
           otherType == SML::SPARSE_MATRIX_ROW);
   
    switch(otherType)
    {
    case SML::DENSE_VECTOR:  output = inner_prod(*sv, *(other.dv)); break;
    case SML::SPARSE_VECTOR: output = inner_prod(*sv, *(other.sv)); break;
    case SML::DENSE_MATRIX_ROW: output = inner_prod(*sv, *(other.dmr)); break;
    case SML::SPARSE_MATRIX_ROW: output = inner_prod(*sv, *(other.smr)); break;
    }  
}


/**   Actual dot product between sv and any dense/sparse matrix.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::sv_dot_m(const TheMatrix& other, TheMatrix& output) const
{
    int otherType  = other.Type();
    int outputType = output.Type();
    
    // type assertion
    assert(otherType == SML::DENSE_MATRIX ||  
           otherType == SML::SPARSE_MATRIX || 
           otherType == SML::DENSE_VECTOR ||
           otherType == SML::SPARSE_VECTOR ||
           otherType == SML::DENSE_MATRIX_ROW ||
           otherType == SML::SPARSE_MATRIX_ROW);
    
    assert(outputType == SML::DENSE_VECTOR ||
           outputType == SML::SPARSE_VECTOR ||
           outputType == SML::DENSE_MATRIX_ROW ||
           outputType == SML::SPARSE_MATRIX_ROW);
    
    switch(otherType) 
    {
    case SML::SPARSE_MATRIX: 
        switch(outputType) 
        {
            //case SML::DENSE_VECTOR:  axpy_prod(*sv, *(other.sm), *(output.dv)); break;
            //case SML::SPARSE_VECTOR: axpy_prod(*sv, *(other.sm), *(output.sv)); break;
            //case SML::DENSE_VECTOR:  axpy_prod(trans(*(other.sm)), *sv, *(output.dv)); break;
        case SML::DENSE_VECTOR:  noalias(*(output.dv)) = prod(*(other.sm), *sv); break;
        case SML::SPARSE_VECTOR: axpy_prod(trans(*(other.sm)), *sv, *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(trans(*(other.sm)), *sv, *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(trans(*(other.sm)), *sv, *(output.smr)); break;
            
        }  
        break;
        
    case SML::DENSE_MATRIX:
        switch(outputType) 
        {
            //case SML::DENSE_VECTOR:  axpy_prod(*sv, *(other.dm), *(output.dv)); break;      
            //case SML::SPARSE_VECTOR: axpy_prod(*sv, *(other.dm), *(output.sv)); break;      
        case SML::DENSE_VECTOR:  axpy_prod(trans(*(other.dm)), *sv, *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(trans(*(other.dm)), *sv, *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(trans(*(other.dm)), *sv, *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(trans(*(other.dm)), *sv, *(output.smr)); break;
        }
        break;
        
        // the following cases are for pathological case: output is a 1x1 matrix
    case SML::DENSE_VECTOR:
    case SML::SPARSE_VECTOR:      
    case SML::DENSE_MATRIX_ROW:
    case SML::SPARSE_MATRIX_ROW:             
        double value = 0.0;
        sv_dot_v(other, value);
        switch(outputType)
        {
        case SML::DENSE_VECTOR: (*output.dv)[0] = value; break;
        case SML::SPARSE_VECTOR: (*output.sv)[0] = value; break;
        case SML::DENSE_MATRIX_ROW: (*output.dmr)[0] = value; break;
        case SML::SPARSE_MATRIX_ROW: (*output.smr)[0] = value; break;                      
        } 
        break;
    }  
}


/**   Actual dot product between smr and any dense/sparse vector.
 *    The result is a double value.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::smr_dot_v(const TheMatrix& other, double& output) const
{
        int otherType = other.Type();
   
        // type assertion
        assert(otherType == SML::DENSE_VECTOR ||
               otherType == SML::SPARSE_VECTOR || 
               otherType == SML::DENSE_MATRIX_ROW || 
               otherType == SML::SPARSE_MATRIX_ROW);
   
        switch(otherType) 
        {
                case SML::DENSE_VECTOR:  output = inner_prod(*smr, *(other.dv)); break;
                case SML::SPARSE_VECTOR: output = inner_prod(*smr, *(other.sv)); break;
                case SML::DENSE_MATRIX_ROW: output = inner_prod(*smr, *(other.dmr)); break;
                case SML::SPARSE_MATRIX_ROW: output = inner_prod(*smr, *(other.smr)); break;
        }  
}


/**   Actual dot product between smr and any dense/sparse matrix.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::smr_dot_m(const TheMatrix& other, TheMatrix& output) const
{
        int otherType  = other.Type();
        int outputType = output.Type();
           
        // type assertion
        assert(otherType == SML::DENSE_MATRIX ||  
               otherType == SML::SPARSE_MATRIX || 
               otherType == SML::DENSE_VECTOR ||
               otherType == SML::SPARSE_VECTOR ||
               otherType == SML::DENSE_MATRIX_ROW ||
               otherType == SML::SPARSE_MATRIX_ROW);
        
        assert(outputType == SML::DENSE_VECTOR ||
               outputType == SML::SPARSE_VECTOR ||
               outputType == SML::DENSE_MATRIX_ROW ||
               outputType == SML::SPARSE_MATRIX_ROW);
   
        switch(otherType)
        {
                case SML::SPARSE_MATRIX: 
                        switch(outputType)
                        {
                                case SML::DENSE_VECTOR: axpy_prod(*smr, *(other.sm), *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(trans(*(other.sm)), *smr, *(output.sv)); break;                                
                                //case SML::SPARSE_VECTOR: noalias(*(output.sv)) = prod(*dv, *(other.sm)); break;
                                case SML::DENSE_MATRIX_ROW: axpy_prod(*smr, *(other.sm), *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(trans(*(other.sm)), *smr, *(output.smr)); break;

                        }  
                        break;
         
                case SML::DENSE_MATRIX:
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*smr, *(other.dm), *(output.dv)); break;      
                                case SML::SPARSE_VECTOR: axpy_prod(*smr, *(other.dm), *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW: axpy_prod(*smr, *(other.dm), *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(trans(*(other.dm)), *smr, *(output.smr)); break;
                        }
                        break;
                
                // the following cases are for pathological case: output is a 1x1 matrix
                case SML::DENSE_VECTOR:
                case SML::SPARSE_VECTOR:      
                case SML::DENSE_MATRIX_ROW:
                case SML::SPARSE_MATRIX_ROW:             
                        double value = 0.0;
                        smr_dot_v(other, value);
                        switch(outputType)
                        {
                                case SML::DENSE_VECTOR: (*output.dv)[0] = value; break;
                                case SML::SPARSE_VECTOR: (*output.sv)[0] = value; break;
                                case SML::DENSE_MATRIX_ROW: (*output.dmr)[0] = value; break;
                                case SML::SPARSE_MATRIX_ROW: (*output.smr)[0] = value; break;                      
                        } 
                        break;
        }
}


/**   Actual dot product between dmr and any dense/sparse vector.
 *    The result is a double value.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::dmr_dot_v(const TheMatrix& other, double& output) const
{
        int otherType = other.Type();
        
        // type assertion
        assert(otherType == SML::DENSE_VECTOR ||
               otherType == SML::SPARSE_VECTOR || 
               otherType == SML::DENSE_MATRIX_ROW || 
               otherType == SML::SPARSE_MATRIX_ROW);
      
        switch(otherType) 
        {
                case SML::DENSE_VECTOR:  output = inner_prod(*dmr, *(other.dv)); break;
                case SML::SPARSE_VECTOR: output = inner_prod(*dmr, *(other.sv)); break;
                case SML::DENSE_MATRIX_ROW: output = inner_prod(*dmr, *(other.dmr)); break;
                case SML::SPARSE_MATRIX_ROW: output = inner_prod(*dmr, *(other.smr)); break;
        }  
}


/**   Actual dot product between dmr and any dense/sparse matrix.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::dmr_dot_m(const TheMatrix& other, TheMatrix& output) const
{
        int otherType  = other.Type();
        int outputType = output.Type();
   
        // type assertion
        assert(otherType == SML::DENSE_MATRIX ||  
               otherType == SML::SPARSE_MATRIX || 
               otherType == SML::DENSE_VECTOR ||
               otherType == SML::SPARSE_VECTOR ||
               otherType == SML::DENSE_MATRIX_ROW ||
               otherType == SML::SPARSE_MATRIX_ROW);
        
        assert(outputType == SML::DENSE_VECTOR ||
               outputType == SML::SPARSE_VECTOR ||
               outputType == SML::DENSE_MATRIX_ROW ||
               outputType == SML::SPARSE_MATRIX_ROW);

   
        switch(otherType)
        {
                case SML::SPARSE_MATRIX: 
                        switch(outputType)
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*dmr, *(other.sm), *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(trans(*(other.sm)), *dmr, *(output.sv)); break;
                                //case SML::SPARSE_VECTOR: noalias(*(output.sv)) = prod(*dv, *(other.sm)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*dmr, *(other.sm), *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(trans(*(other.sm)), *dmr, *(output.smr)); break;
                        }  
                        break;
         
                case SML::DENSE_MATRIX:
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*dmr, *(other.dm), *(output.dv)); break;      
                                case SML::SPARSE_VECTOR: axpy_prod(*dmr, *(other.dm), *(output.sv)); break;      
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*dmr, *(other.dm), *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(trans(*(other.dm)), *dmr, *(output.smr)); break;
                        }
                        break;
                        
                // the following cases are for pathological case: output is a 1x1 matrix
                case SML::DENSE_VECTOR:
                case SML::SPARSE_VECTOR:      
                case SML::DENSE_MATRIX_ROW:
                case SML::SPARSE_MATRIX_ROW:             
                        double value = 0.0;
                        dmr_dot_v(other, value);
                        switch(outputType)
                        {
                                case SML::DENSE_VECTOR: (*output.dv)[0] = value; break;
                                case SML::SPARSE_VECTOR: (*output.sv)[0] = value; break;
                                case SML::DENSE_MATRIX_ROW: (*output.dmr)[0] = value; break;
                                case SML::SPARSE_MATRIX_ROW: (*output.smr)[0] = value; break;                      
                        } 
                        break;
   }
}


/**   Actual dot product between dm and any dense/sparse matrix/vector.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::dm_dot(const TheMatrix& other, TheMatrix& output) const
{
    int otherType  = other.Type();
    int outputType = output.Type();
    
    // type assertion
    assert(otherType == SML::DENSE_MATRIX ||  
           otherType == SML::SPARSE_MATRIX || 
           otherType == SML::DENSE_VECTOR ||
           otherType == SML::SPARSE_VECTOR ||
           otherType == SML::DENSE_MATRIX_ROW ||
           otherType == SML::SPARSE_MATRIX_ROW);
    
    assert(outputType == SML::DENSE_VECTOR ||
           outputType == SML::SPARSE_VECTOR ||
           outputType == SML::DENSE_MATRIX_ROW ||
           outputType == SML::SPARSE_MATRIX_ROW);
    
    switch(otherType) 
    {
    case SML::DENSE_VECTOR: 
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
        switch(outputType) 
        {
        case SML::DENSE_VECTOR:  axpy_prod(*dm, *(other.dv), *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(*dm, *(other.dv), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*dm, *(other.dv), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*dm, *(other.dv), *(output.smr)); break;
        }
        break;
        
    case SML::SPARSE_VECTOR: 
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
        switch(outputType) 
        {
        case SML::DENSE_VECTOR:  axpy_prod(*dm, *(other.sv), *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(*dm, *(other.sv), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*dm, *(other.sv), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*dm, *(other.sv), *(output.smr)); break;    
        }
        break;
        
    case SML::DENSE_MATRIX: 
        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
        switch(outputType)
        {
        case SML::DENSE_MATRIX:  axpy_prod(*dm, *(other.dm), *(output.dm)); break;
        case SML::SPARSE_MATRIX: noalias(*(output.sm)) = prod(*dm, *(other.dm)); break;
            //axpy_prod(*dm, *(other.dm), *(output.sm)); break;
            // prod is faster in this case	       
        }
        break;
        
    case SML::SPARSE_MATRIX: 
        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
        switch(outputType) 
        {
        case SML::DENSE_MATRIX:  axpy_prod(*dm, *(other.sm), *(output.dm)); break;
        case SML::SPARSE_MATRIX: axpy_prod(*dm, *(other.sm), *(output.sm)); break;	       
        }
        break;
        
    case SML::DENSE_MATRIX_ROW:
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);                        
        switch(outputType) 
        {
        case SML::DENSE_VECTOR:  axpy_prod(*dm, *(other.dmr), *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(*dm, *(other.dmr), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*dm, *(other.dmr), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*dm, *(other.dmr), *(output.smr)); break;                               
        }
        break;
        
    case SML::SPARSE_MATRIX_ROW:
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
        switch(outputType) 
        {
        case SML::DENSE_VECTOR: axpy_prod(*dm, *(other.smr), *(output.dv)); break;
            //noalias(*(output.dv)) = prod(*dm, *(other.smr)); break; // this is slower!
        case SML::SPARSE_VECTOR: axpy_prod(*dm, *(other.smr), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*dm, *(other.smr), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*dm, *(other.smr), *(output.smr)); break;                                
        }
        break;
    }  
}


/**   Actual dot product between sm and any dense/sparse matrix/vector.
 * 
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::sm_dot(const TheMatrix& other, TheMatrix& output) const
{
    int otherType  = other.Type();
    int outputType = output.Type();
    
    // type assertion
    assert(otherType == SML::DENSE_MATRIX ||  
           otherType == SML::SPARSE_MATRIX || 
           otherType == SML::DENSE_VECTOR ||
           otherType == SML::SPARSE_VECTOR ||
           otherType == SML::DENSE_MATRIX_ROW ||
           otherType == SML::SPARSE_MATRIX_ROW);
    
    assert(outputType == SML::DENSE_VECTOR ||
           outputType == SML::SPARSE_VECTOR ||
           outputType == SML::DENSE_MATRIX_ROW ||
           outputType == SML::SPARSE_MATRIX_ROW);
    
    switch(otherType) 
    {
    case SML::DENSE_VECTOR: 
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
        switch(outputType) 
        {
        case SML::DENSE_VECTOR:  axpy_prod(*sm, *(other.dv), *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(*sm, *(other.dv), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*sm, *(other.dv), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*sm, *(other.dv), *(output.smr)); break;
        }
        break;
        
    case SML::SPARSE_VECTOR: 
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
        switch(outputType) 
        {
        case SML::DENSE_VECTOR:  axpy_prod(*sm, *(other.sv), *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(*sm, *(other.sv), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*sm, *(other.sv), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*sm, *(other.sv), *(output.smr)); break;
        }
        break;
        
    case SML::DENSE_MATRIX: 
        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
        switch(outputType) 
        {
        case SML::DENSE_MATRIX:  axpy_prod(*sm, *(other.dm), *(output.dm)); break;
        case SML::SPARSE_MATRIX: noalias(*(output.sm)) = prod(*(sm), *(other.dm)); break;
            // axpy_prod(*sm, *(other.dm), *(output.sm)); break;
            // prod is about 2 times faster in this case
        }
        break;
        
        // TODO: try block_prod for sparse matrix multiplication
    case SML::SPARSE_MATRIX: 
        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
        switch(outputType) 
        {
        case SML::DENSE_MATRIX: axpy_prod(*sm, *(other.sm), *(output.dm)); break;
            //noalias(*(output.dm)) = prod(*sm, *(other.sm)); break;
        case SML::SPARSE_MATRIX: noalias(*(output.sm)) = prod(*(sm), *(other.sm)); break;
            //axpy_prod(*sm, *(other.sm), *(output.sm)); break;
            // prod is slightly faster than axpy
        }
        break;
        
    case SML::DENSE_MATRIX_ROW:
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
        switch(outputType) 
        {
        case SML::DENSE_VECTOR:  axpy_prod(*sm, *(other.dmr), *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(*sm, *(other.dmr), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*sm, *(other.dmr), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*sm, *(other.dmr), *(output.smr)); break;
        }
        break;
        
    case SML::SPARSE_MATRIX_ROW:
        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
        switch(outputType) 
        {
        case SML::DENSE_VECTOR:  axpy_prod(*sm, *(other.smr), *(output.dv)); break;
        case SML::SPARSE_VECTOR: axpy_prod(*sm, *(other.smr), *(output.sv)); break;
        case SML::DENSE_MATRIX_ROW: axpy_prod(*sm, *(other.smr), *(output.dmr)); break;
        case SML::SPARSE_MATRIX_ROW: axpy_prod(*sm, *(other.smr), *(output.smr)); break;
        }
        break;
    }    
}


/**   Actual dot product between transposed sm and any dense/sparse matrix/vector.
 *    (uBlas will do the matrix size checking in debug mode)
 *
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::sm_transposedot(const TheMatrix& other, TheMatrix& output) const
{
#if DEBUG
        printf("In sm_transposedot(const TheMatrix&, TheMatrix&) const");
#endif  
    
        int otherType  = other.Type();
        int outputType = output.Type();
   
        // type assertion
        assert(otherType == SML::DENSE_MATRIX ||  
               otherType == SML::SPARSE_MATRIX || 
               otherType == SML::DENSE_VECTOR ||
               otherType == SML::SPARSE_VECTOR ||
               otherType == SML::DENSE_MATRIX_ROW ||
               otherType == SML::SPARSE_MATRIX_ROW);
        
        assert(outputType == SML::DENSE_VECTOR ||
               outputType == SML::SPARSE_VECTOR ||
               outputType == SML::DENSE_MATRIX_ROW ||
               outputType == SML::SPARSE_MATRIX_ROW);
        
        switch(otherType) 
        {
                case SML::DENSE_VECTOR:    
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.dv), *sm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.dv), *sm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.dv), *sm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.dv), *sm, *(output.smr)); break;
                        }
                        break;
	 
                case SML::SPARSE_VECTOR: 
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.sv), *sm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.sv), *sm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.sv), *sm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.sv), *sm, *(output.smr)); break;
                        }
                        break;
	 
                case SML::DENSE_MATRIX: 
                        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
                        switch(outputType) 
                        {
                                case SML::DENSE_MATRIX:  axpy_prod(ublas::trans(*sm), *(other.dm), *(output.dm)); break;
                                case SML::SPARSE_MATRIX: noalias(*(output.sm)) = prod(ublas::trans(*sm), *(other.dm)); break;
                        }
                        break;
	 
                case SML::SPARSE_MATRIX: 
                        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
                        switch(outputType) 
                        {
                                case SML::DENSE_MATRIX:  axpy_prod(ublas::trans(*sm), *(other.sm), *(output.dm)); break;                 
                                case SML::SPARSE_MATRIX: noalias(*(output.sm)) = prod(ublas::trans(*sm), *(other.sm)); break;
                        }
                        break;

                case SML::DENSE_MATRIX_ROW:
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                        switch(outputType) 
                        {                                
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.dmr), *sm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.dmr), *sm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.dmr), *sm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.dmr), *sm, *(output.smr)); break;
                        }
                        break;

                case SML::SPARSE_MATRIX_ROW:
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.smr), *sm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.smr), *sm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.smr), *sm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.smr), *sm, *(output.smr)); break;
                        }
                        break;
        }  
}



/**   Actual dot product between transposed dm and any dense/sparse matrix/vector.
 *    (uBlas will do the matrix size checking in debug mode)
 *
 *    \param other  - [read]  The left operand of the dot product
 *    \param output - [write] The result of the dot product
 */
void TheMatrix::dm_transposedot(const TheMatrix& other, TheMatrix& output) const
{
#if DEBUG
        printf("In sm_transposedot(const TheMatrix&, TheMatrix&) const");
#endif  

        int otherType  = other.Type();
        int outputType = output.Type();
   
        // type assertion
        assert(otherType == SML::DENSE_MATRIX ||  
               otherType == SML::SPARSE_MATRIX || 
               otherType == SML::DENSE_VECTOR ||
               otherType == SML::SPARSE_VECTOR ||
               otherType == SML::DENSE_MATRIX_ROW ||
               otherType == SML::SPARSE_MATRIX_ROW);
        
        assert(outputType == SML::DENSE_VECTOR ||
               outputType == SML::SPARSE_VECTOR ||
               outputType == SML::DENSE_MATRIX_ROW ||
               outputType == SML::SPARSE_MATRIX_ROW);
        
        switch(otherType) 
        {
                case SML::DENSE_VECTOR: 
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.dv), *dm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.dv), *dm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.dv), *dm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.dv), *dm, *(output.smr)); break;
                        }
                        break;
	 
                case SML::SPARSE_VECTOR: 
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                	switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.sv), *dm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.sv), *dm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.sv), *dm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.sv), *dm, *(output.smr)); break;
                        }
                        break;
	 
                case SML::DENSE_MATRIX: 
                        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
                        switch(outputType) 
                        {
                                case SML::DENSE_MATRIX:  axpy_prod(ublas::trans(*dm), *(other.dm), *(output.dm)); break;
                                case SML::SPARSE_MATRIX: noalias(*(output.sm)) = prod(ublas::trans(*dm), *(other.dm)); break;
                        }
                        break;
	 
                case SML::SPARSE_MATRIX: 
                        assert(outputType == SML::DENSE_MATRIX  ||  outputType == SML::SPARSE_MATRIX);
                        switch(outputType) 
                        {
                                case SML::DENSE_MATRIX:  axpy_prod(ublas::trans(*dm), *(other.sm), *(output.dm)); break;                 
                                case SML::SPARSE_MATRIX: noalias(*(output.sm)) = prod(ublas::trans(*dm), *(other.sm)); break;
                        }
                        break;         

                case SML::DENSE_MATRIX_ROW:
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.dmr), *dm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.dmr), *dm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.dmr), *dm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.dmr), *dm, *(output.smr)); break;
                        }
                        break;

                case SML::SPARSE_MATRIX_ROW:
                        assert(outputType == SML::DENSE_VECTOR  ||  outputType == SML::SPARSE_VECTOR || outputType == SML::DENSE_MATRIX_ROW || outputType == SML::SPARSE_MATRIX_ROW);
                        switch(outputType) 
                        {
                                case SML::DENSE_VECTOR:  axpy_prod(*(other.smr), *dm, *(output.dv)); break;
                                case SML::SPARSE_VECTOR: axpy_prod(*(other.smr), *dm, *(output.sv)); break;
                                case SML::DENSE_MATRIX_ROW:  axpy_prod(*(other.smr), *dm, *(output.dmr)); break;
                                case SML::SPARSE_MATRIX_ROW: axpy_prod(*(other.smr), *dm, *(output.smr)); break;
                        }
                        break;
        }
   
}


/**   Actual function to add/set a dense row of values into sparse/dense vector. 
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing), must be 0
 *    \param length  - [read] The number of elements in that row
 *    \param values  - [read] The row elements 
 */
void TheMatrix::v_addrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values)
{
        assert(type == SML::DENSE_VECTOR || type == SML::SPARSE_VECTOR);
        assert(rowidx == 0);
        assert(length == this->length);
        assert(values);
        
        switch(type)
        {
                case SML::DENSE_VECTOR:
                        for(unsigned int i=0; i<length; i++)
                                (*dv)[i] = values[i];
                        break;
                case SML::SPARSE_VECTOR:
                        for(unsigned int i=0; i<length; i++)
                                (*sv)[i] = values[i];
                        break;
        }   
}

/**   Actual function to add/set a sparse row of values into sparse/dense vector. 
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing), must be 0
 *    \param length  - [read] The number of nnz elements in that row
 *    \param values  - [read] The nonzero elements
 *    \param indices - [read] The column indices corresponding to the nonzero elements
 */
void TheMatrix::v_addrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices)
{
        assert(type == SML::DENSE_VECTOR || type == SML::SPARSE_VECTOR);
        assert(rowidx == 0);
        assert(length <= this->length);
        assert(values);
        assert(indices);
        
        switch(type)
        {
                case SML::DENSE_VECTOR:
                        for(unsigned int i=0; i<length; i++)
                                (*dv)[indices[i]] = values[i];
                        break;
                case SML::SPARSE_VECTOR:
                        for(unsigned int i=0; i<length; i++)
                                (*sv)[indices[i]] = values[i];
                        break;
        }
}


/**   Actual function to get a dense row of values from dense/sparse vector.
 *
 *    \param rowidx  - [read]  The row which returning values belong to (zero-based indexing)
 *    \param length  - [write] The nnz of the row
 *    \param values  - [write] The values in the row (memory preallocated) 
 */
void TheMatrix::v_getrow_dense(const unsigned int& rowidx, unsigned int& length, double* values) const
{
        assert(type == SML::DENSE_VECTOR || type == SML::SPARSE_VECTOR);
        assert(rowidx == 0);
        assert(values);
        
        switch(type)
        {
                case SML::DENSE_VECTOR: 
                        for(unsigned int i=0; i<this->length; i++)
                                values[i] = (*dv)[i];
                        break;
                case SML::SPARSE_VECTOR: 
                        for(unsigned int i=0; i<this->length; i++)
                                values[i] = (*sv)[i];
                        break;                      
        }   
        length = this->length;
}



/**   Actual function to get a sparse row of values from dense/sparse vector.
 *    This function is provide for consistency purpose.
 *
 *    \param rowidx  - [read]  The row which returning values belong to (zero-based indexing)
 *    \param length  - [write] The nnz of the row
 *    \param values  - [write] The values in the row (memory preallocated)
 *    \param indices - [write] The column indices of the values in the row (memory preallocated)
 */
void TheMatrix::v_getrow_sparse(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int* indices) const
{
        assert(type == SML::DENSE_VECTOR || type == SML::SPARSE_VECTOR);
        assert(rowidx == 0);
        assert(values);
        assert(indices);

        int nnz = 0;
        switch(type)
        {
                case SML::DENSE_VECTOR: 
                        for(unsigned int i=0; i<this->length; i++)
                                if(SML::abs((*dv)[i]) > 0.0)
                                {
                                        values[nnz] = (*dv)[i];        
                                        indices[nnz] = i;
                                        nnz++;
                                }
                                
                        break;
                case SML::SPARSE_VECTOR: 
                        for(unsigned int i=0; i<this->length; i++)
                                //if(SML::abs((*sv)[i]) > 0.0)  // not supported?
                                if((*sv)[i] > 0.0 || (*sv)[i] < 0.0)
                                {
                                        values[nnz] = (*sv)[i];        
                                        indices[nnz] = i;
                                        nnz++;
                                }                                
                        break;                      
        }   
        length = nnz;
}



/**   Actual function to add/set a dense row of values into dense matrix. 
 *    Accept sparse or dense array representation.
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing)
 *    \param length  - [read] The number of elements in that row
 *    \param values  - [read] The row elements
 */
void TheMatrix::dm_addrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values)
{
        assert(type == SML::DENSE_MATRIX);
        assert(length == this->column);
        assert(values);

        for(unsigned int i=0; i<length; i++)        
                (*dm)(rowidx,i) = values[i];      
}


/**   Actual function to add/set a sparse row of values into dense matrix. 
 *    Accept sparse or dense array representation.
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing)
 *    \param length  - [read] The number of elements in that row
 *    \param values  - [read] The row elements
 *    \param indices - [read] The column indices corresponding to the values 
 */
void TheMatrix::dm_addrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices)
{
        assert(type == SML::DENSE_MATRIX);
        assert(length <= this->column);
        assert(values);
        assert(indices);

        for(unsigned int i=0; i<length; i++)
                (*dm)(rowidx,indices[i]) = values[i];

}


/**   Actual function to add a dense row of values into sparse matrix.
 *    Only for (sequential) matrix initialization, i.e. m(i,j) = value \forall i \in [row] and j \in [col]
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing)
 *    \param length  - [read] The number of elements in that row
 *    \param values  - [read] The row elements
 */
void TheMatrix::sm_addrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values)
{
        assert(type == SML::SPARSE_MATRIX);   
        assert(length == this->column);   
        assert(values);

        for(unsigned int i=0; i<length; i++)
                if(SML::abs(values[i]) > ZERO_EPS)
                        sm->push_back(rowidx, i, values[i]);
   
}


/**   Actual function to add a sparse row of values into sparse matrix.
 *    Only for (sequential) matrix initialization, i.e. m(i,j) = value \forall i \in [row] and j \in [col]
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing)
 *    \param length  - [read] The number of elements in that row
 *    \param values  - [read] The row elements
 *    \param indices - [read] The column indices corresponding to the values 
 */
void TheMatrix::sm_addrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices)
{
        assert(type == SML::SPARSE_MATRIX);
        assert(length <= this->column);
        assert(values);
        assert(indices);

        for(unsigned int i=0; i<length; i++)
                sm->push_back(rowidx, indices[i], values[i]);   
}


/**   Actual function to set a dense row of values into sparse matrix.
 *    Can be used at anytime! (but slower than sm_addrow when used to initialize a matrix in sequential order)
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing)
 *    \param length  - [read] The number of elements in that row
 *    \param values  - [read] The row elements
 */
void TheMatrix::sm_setrow_dense(const unsigned int& rowidx, const unsigned int& length, const double* values)
{
        assert(type == SML::SPARSE_MATRIX);
        assert(length == this->column);
        assert(values);

        for(unsigned int i=0; i<length; i++) 
                if(SML::abs(values[i]) > ZERO_EPS)
                        (*sm)(rowidx, i) = values[i];
}


/**   Actual function to set a sparse row of values into sparse matrix.
 *    Can be used at anytime! (but slower than sm_addrow when used to initialize a matrix in sequential order)
 *
 *    \param rowidx  - [read] The row which incoming values belong to (zero-based indexing)
 *    \param length  - [read] The number of elements in that row
 *    \param values  - [read] The row elements
 *    \param indices - [read] The column indices corresponding to the values 
 */
void TheMatrix::sm_setrow_sparse(const unsigned int& rowidx, const unsigned int& length, const double* values, const unsigned int* indices)
{
        assert(type == SML::SPARSE_MATRIX);
        assert(length <= this->column);
        assert(values);
        assert(indices);

        for(unsigned int i=0; i<length; i++) 
                (*sm)(rowidx, indices[i]) = values[i];
}


/**   
 * Actual function to get a dense row of values from dense matrix. NOTE:
 * The values array is allocated by this function.  It is the duty of
 * the caller to free this memory.
 * 
 * CAUTION: This function works on the internal structures of uBLAS
 *
 *    \param rowidx  - [read]  The row which returning values belong to (zero-based indexing)
 *    \param length  - [write] The number of elements in the row
 *    \param values  - [write] The values in the row (memory preallocated)
 */
void TheMatrix::dm_getrow_dense(const unsigned int& rowidx, unsigned int& length, double* values) const
{
        assert(type == SML::DENSE_MATRIX);
        assert(values);
        assert(rowidx < row);
  
        // Assumption:
        // 1. uBLAS matrix store matrix elements in a 1D array
        // 2. matrix element is of type read (as #defined)
  
        unsigned int start_idx = rowidx*column;
        length = this->column;        
        memcpy(values, &((*dm).data()[start_idx]), sizeof(double)*length);        
}


/**   Actual function to get a sparse row of values from dense matrix.
 *    (Note that this function works on the internal structures of uBLAS matrix)
 *
 *    \param rowidx  - [read]  The row which returning values belong to (zero-based indexing)
 *    \param length  - [write] The nnz of the row
 *    \param values  - [write] The values in the row (memory preallocated)
 *    \param indices - [write] The column indices of the values in the row (memory preallocated)
 */
void TheMatrix::dm_getrow_sparse(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int* indices) const
{
        assert(type == SML::DENSE_MATRIX);
        assert(values);
        assert(indices);
        assert(rowidx < row);

        // Assumption:
        // 1. uBLAS matrix store matrix elements in a 1D array

        unsigned int start_idx = rowidx*this->column;
        double tmp_value = 0;
        length = 0;

        for(unsigned int i=0; i<this->column; i++) 
        {
                tmp_value = (*dm).data()[start_idx+i];
                if(SML::abs(tmp_value) > ZERO_EPS) 
                {
                        values[length]  = tmp_value;
                        indices[length] = i;
                        length++;
                }
        }
}



/**   Actual function to get a dense row of values from sparse matrix.
 *    (Note that this function works on the internal structures of uBLAS compressed_matrix)
 *
 *    \param rowidx  - [read]  The row which returning values belong to (zero-based indexing)
 *    \param length  - [write] The number of elements in the row
 *    \param values  - [write] The values in the row (memory preallocated)
 */
void TheMatrix::sm_getrow_dense(const unsigned int& rowidx, unsigned int& length, double* values) const
{
        assert(type == SML::SPARSE_MATRIX);
        assert(values);
        assert(rowidx < row);

        // Assumption:
        // 1. uBLAS compressed_matrix is of standard (netlib.org) compressed sparse row format

        length = this->column;
        memset(values,0,sizeof(double)*length);

        unsigned int row_nnz = (*sm).index1_data()[rowidx+1] - (*sm).index1_data()[rowidx];
        unsigned int start_idx = (*sm).index1_data()[rowidx];
        unsigned int index = 0;
        double tmp_value = 0;

        for(unsigned int i=0; i<row_nnz; i++) 
        {
                index = (*sm).index2_data()[start_idx+i];
                tmp_value = (*sm).value_data()[start_idx+i];
                values[index] = tmp_value;
        }
}


/**   Actual function to get a sparse row of values from sparse matrix.
 *    (Note that this function works on the internal structures of uBLAS compressed_matrix)
 *
 *    \param rowidx  - [read]  The row which returning values belong to (zero-based indexing)
 *    \param length  - [write] The nnz of the row
 *    \param values  - [write] The values in the row (memory preallocated)
 *    \param indices - [write] The column indices of the values in the row (memory preallocated)
 */
void TheMatrix::sm_getrow_sparse(const unsigned int& rowidx, unsigned int& length, double* values, unsigned int* indices) const
{
        assert(type == SML::SPARSE_MATRIX);
        assert(values);
        assert(indices);
        assert(rowidx < row);

        // Assumption:
        // 1. uBLAS compressed_matrix is of standard (netlib.org) compressed sparse row format
 
        length = (*sm).index1_data()[rowidx+1] - (*sm).index1_data()[rowidx];
        unsigned int start_idx = (*sm).index1_data()[rowidx];
        for(unsigned int i=0; i<length; i++) 
        {
                values[i] = (*sm).value_data()[start_idx+i];
                indices[i] = (*sm).index2_data()[start_idx+i];
        }
}


/**   Create a new TheMatrix referencing to the index-th row of *this.
 *
 *    \param index [read] The row number
 */
TheMatrix* TheMatrix::CreateMatrixRowView(const unsigned int& index)
{
        assert(type == SML::SPARSE_MATRIX || type == SML::DENSE_MATRIX);
   
        TheMatrix *mat_row = new TheMatrix();
        mat_row->isReference = true;

        switch(type) 
        {
                case SML::SPARSE_MATRIX:
                        mat_row->smr = new ublas::matrix_row<ublas::compressed_matrix<double> >(*sm,index);
                        mat_row->type = SML::SPARSE_MATRIX_ROW;
                        mat_row->row = 1;
                        mat_row->column = mat_row->smr->size();
                        mat_row->dense_or_sparse = SML::SPARSE;
                        break;
         
                case SML::DENSE_MATRIX:
                        mat_row->dmr = new ublas::matrix_row<ublas::matrix<double> >(*dm,index);
                        mat_row->type = SML::DENSE_MATRIX_ROW;
                        mat_row->row = 1;
                        mat_row->column = mat_row->dmr->size();
                        mat_row->dense_or_sparse = SML::DENSE;
                        break;
        }
        mat_row->vector_or_matrix = SML::VECTOR;
        mat_row->column_or_row = SML::ROW;
        mat_row->assign_pointers();

        return mat_row;
}

/**   Assign a reference to the index-th row of *this to mat_row
 *
 *    @param index [read] The row number
 *    @param mat_row [write] The matrix row view
 */
void TheMatrix::CreateMatrixRowView(const unsigned int& index, TheMatrix* mat_row)
{
        assert(type == SML::SPARSE_MATRIX || type == SML::DENSE_MATRIX);
		assert(mat_row);
		
        //TheMatrix *mat_row = new TheMatrix();
        mat_row->isReference = true;

        switch(type) 
        {
                case SML::SPARSE_MATRIX:
					if(mat_row->smr)
						delete mat_row->smr;
					mat_row->smr = new ublas::matrix_row<ublas::compressed_matrix<double> >(*sm,index);
					mat_row->type = SML::SPARSE_MATRIX_ROW;
					mat_row->row = 1;
					mat_row->column = mat_row->smr->size();
					mat_row->dense_or_sparse = SML::SPARSE;
					break;
					
                case SML::DENSE_MATRIX:
					if(mat_row->dmr)
						delete mat_row->dmr;
					mat_row->dmr = new ublas::matrix_row<ublas::matrix<double> >(*dm,index);
					mat_row->type = SML::DENSE_MATRIX_ROW;
					mat_row->row = 1;
					mat_row->column = mat_row->dmr->size();
					mat_row->dense_or_sparse = SML::DENSE;
					break;
        }
        mat_row->vector_or_matrix = SML::VECTOR;
        mat_row->column_or_row = SML::ROW;
        mat_row->assign_pointers();
}

#endif
