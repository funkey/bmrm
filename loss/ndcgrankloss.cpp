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
 *
 * Created: (06/01/2008) 
 *
 * Last Updated:
 */

#ifndef _NDCGRANKLOSS_CPP_
#define _NDCGRANKLOSS_CPP_

#include "ndcgrankloss.hpp"
#include "sml.hpp"
#include "lap.hpp"
#include <ext/numeric>         // for iota
#include "configuration.hpp"

using namespace std;

CNDCGRankLoss::CNDCGRankLoss(CModel* &model, CVecData* &data)
   : CRankLoss(model, data) 
{
   Configuration &config = Configuration::GetInstance();
   scalingFactor = 1.0/(10.0*m);
   
   if(config.IsSet("NDCG.truncation")) 
      truncation = config.GetInt("NDCG.truncation"); 
   else
      truncation = 10; /* default value is 10 */
   
  double c_value = 0;
  if(config.IsSet("NDCG.c_function")) 
     c_value = config.GetDouble("NDCG.c_function"); 
  
  cout << "In CNDCGRankLoss, truncation = "<< truncation << ", "<<"c_value = "<< c_value <<endl;
  int num_query = _data->NumOfSubset();
  double *Y_array = _data->labels().Data();
  max_subset_size = 0;
  
  /* sort the documents in each query */
  for (int q=0;q<num_query;q++)
  {
     int offset = _data->subset[q].startIndex;
     int subsetsize = _data->subset[q].size;
     if (max_subset_size < subsetsize)
        max_subset_size = subsetsize;

     int length = subsetsize;
     vector<int> indices(length);
     iota(indices.begin(), indices.end(), 0);
     indirect_greater_than<double> igt(&(Y_array[offset]));
     sort(indices.begin(), indices.end(),igt);
     sort_vectors.push_back(indices);
     /* use this to compute bs */
     vector<double> b;
     compute_b(offset, subsetsize, Y_array, indices, b);
     bs.push_back( b );
  }

  /* use largest query information to compute c_i */
  for(int i=0;i<max_subset_size;i++)
  {
     if(i < truncation)
     {
        if (c_value == 0)
           c.push_back( max(1.0*log(2)/log(i+2), 0.0));
        else if (c_value == 100)
           c.push_back(max(truncation + 1 - i,0));
        else
           c.push_back( max(1.0/pow((i+1),c_value), 0.0));
     }
     else
     {
        c.push_back(0);
     }
  }
  pi = (int *) malloc(sizeof(int)*max_subset_size);
  compute_a(max_subset_size, a); 
}



void CNDCGRankLoss::compute_b(int offset, int size, double *y, 
                              vector<int> ideal_pi, 
                              vector<double> &input_b)
{
   double dy;
   eval_dcg(offset, size, truncation, y, ideal_pi, dy);
   for (int i=0;i<size;i++)
   {
      if (dy == 0)
         input_b.push_back(1);
      else
         input_b.push_back( (pow(2.0,(double) y[offset + ideal_pi[i]]) - 1.0)/dy );    
   }
}

void CNDCGRankLoss::compute_a(int size, vector<double> &input_a)
{
   for (int i=0;i<size;i++)
      if (i<truncation)
         input_a.push_back(log(2)/log(i+2));
      else
         input_a.push_back(0.0);
}

void CNDCGRankLoss::find_permutation(int size, int offset, vector<double> a, 
                                     vector<double> b, vector<double> c, 
                                     double *f, int *input_pi)
{
   /* setting up C_ij */
   int minsize = min(truncation,size);
   double **C = (double **)malloc(sizeof(double*)*minsize);
   for (int i=0;i<minsize;i++)
      C[i] = (double *)malloc(sizeof(double)*size);
   for (int i=0;i<minsize;i++)
      for (int j=0;j<size;j++)
         C[i][j] = -c[i]*get(f, offset, j) + b[j]*a[i];
   
   int *row = (int *)malloc(sizeof(int)*size);
   
   lap(size, minsize, C, input_pi, row);  
   //lap(size, minsize, C, input_pi, row);  
   
   free(row);
   for (int i=0;i<minsize;i++)
      free(C[i]);
   free(C);
}

void CNDCGRankLoss::delta(int size, vector<double> a, vector<double> b, 
			  int *pi, double &value)
{
   value = 0;
   for (int i=0;i<size;i++)
      value += (b[i]- b[pi[i]])*a[i];
}


/**  
 *  Compute NDCGRank loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CNDCGRankLoss::Loss(double& loss, TheMatrix& f)
{
   loss = 0.0;	
   double* f_array = f.Data();  
   for(unsigned int q=0; q < _data->NumOfSubset(); q++)
   {
      int offset = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      current_ideal_pi = sort_vectors[q];
      vector<double> b = bs[q];
      int size  = subsetsize > truncation? truncation: subsetsize;
      
      /* find the best permutation */
      find_permutation(size, offset, a, b, c, f_array, pi);
      
      /* compute the loss */
      double value;
      delta(size, a, b, pi, value);
      
      loss += value;
      
      for (int i=0;i<size;i++)
         loss = loss + c[i]*(get(f_array, offset, pi[i]) - get(f_array, offset, i));
    }
}


/**  
 *  Compute loss and partial derivative of NDCGRank loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CNDCGRankLoss::LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l)
{
   TheMatrix predict;
   loss = 0.0;	
   l.Zero();  
   double* f_array = f.Data();  
   for(unsigned int q=0; q < _data->NumOfSubset(); q++)
   {
      int offset = _data->subset[q].startIndex;
      int subsetsize = _data->subset[q].size;
      current_ideal_pi = sort_vectors[q];
      vector<double> b = bs[q];
      int size  = subsetsize > truncation? truncation: subsetsize;
      //int size = subsetsize;
      
      /* find the best permutation */
      find_permutation(size, offset, a, b, c, f_array, pi);

      /* compute the loss */
      double value;
      //delta(subsetsize, a, b, pi, value);
      delta(size, a, b, pi, value);

      loss += value;
      
      //for (int i=0;i<subsetsize;i++){
      for (int i=0;i<size;i++)
         loss = loss + c[i]*(get(f_array, offset, pi[i]) - get(f_array, offset, i));

      for (int i=0;i<size;i++)
      {
         add(l, offset, i, - c[i]);
         add(l, offset, pi[i], c[i]);
      }
   }
}


double CNDCGRankLoss::get(double *f, int offset, int i)
{
   return f[offset + current_ideal_pi[i]];
   //return f[offset + i];
}

void CNDCGRankLoss::add(TheMatrix &l, int offset, int i, double value)
{
   double temp;
   l.Get(offset + current_ideal_pi[i], temp);
   //l.Get(offset + i, temp);
   l.Set(offset + current_ideal_pi[i], temp + value);
   //l.Set(offset + i, temp + value);
}

#endif
