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
 * Authors: Quoc Viet Le (quocle@stanford.edu)
 *          Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (02/11/2007) 
 *
 * Last Updated: (13/11/2007)   
 */

#include <cmath>
#include "rankloss.hpp"
#include "configuration.hpp"

using namespace std;


/** Transform function value f := <w,x> into label
 *
 *  @param f [read/write] function values / labels
 */
void CRankLoss::Transform(TheMatrix& f)
{
   double* f_array = f.Data(); 
   int len = f.Length();
   for(int i=0; i < len; i++)
      f_array[i] = SML::sgn(f_array[i]);
}

void CRankLoss::eval_dcg(int offset, int size, int truncation, 
                         double *y, vector<int> pi, double &dcg)
{
   int k = size > truncation? truncation: size;
   dcg = 0;
   for (int i=0;i<k;i++)
      dcg += log(2.0)*(pow(2.0, (double)y[offset + pi[i]])-1)/log(i+2.0);
}

void CRankLoss::eval_ndcg(int offset, int size, int truncation, 
                          double *y, double *f, double &ndcg)
{
   vector<int> ideal_pi(size);
   for(size_t i=0; i<ideal_pi.size(); i++) ideal_pi[i] = i;
   
   indirect_greater_than<double> igt1(&(y[offset]));
   sort(ideal_pi.begin(), ideal_pi.end(), igt1);  // ascending order
   
   vector<int> pi(size);
   for(size_t i=0; i<pi.size(); i++) pi[i] = i;
   
   indirect_greater_than<double> igt2(&(f[offset]));
   sort(pi.begin(), pi.end(), igt2);  // ascending order
   
   double ideal_dcg;
   eval_dcg(offset, size, truncation, y, ideal_pi, ideal_dcg);
   
   double dcg;
   eval_dcg(offset, size, truncation, y, pi, dcg);
   
   if (ideal_dcg != 0)
      ndcg = dcg / ideal_dcg;
   else
      ndcg = 1; // ndcg is 1 if the query contains all documents rated 0
}


void CRankLoss::eval_auc(int offset, int size, int truncation, 
                         double *y, double *f, double &auc)
{
   // START : auc
   auc = 0;
   vector<int> idx(size);
   for(size_t i=0; i<idx.size(); i++) idx[i] = i;
   
   indirect_less_than<double> ilt(&(f[offset]));
   sort(idx.begin(), idx.end(), ilt);  // ascending order
   
   // count swapped pairs i.e., misclassified points
   long n_pos = 0, n_neg = 0;
   double sum = 0.0;
   for(int i = 0; i < size; i++) 
   {
      if(y[offset + idx[i]] < 0) 
      {
         sum += n_pos;
         n_neg += 1;
      }
      else
         n_pos += 1;
   }
   if(n_pos == 0 || n_neg == 0)
      auc = 0.0;
   else
      auc = 1.0 - sum/n_neg/n_pos;
}



/** Evaluate the performance of the model on _data
 *
 *  @param f [read] function values
 *  @param predict [write] prdicted labels
 */
void CRankLoss::Perf(TheMatrix& f, TheMatrix& predict)
{
   double* y_array = _data->labels().Data();
   double* f_array = f.Data();
   int max_truncation = 10;
   
   cout << "Performance:" << endl;;
   for(int truncation=1; truncation <= max_truncation; truncation++)
   {
      double mean_ndcg = 0;
      double mean_auc = 0;
      int num_query = _data->NumOfSubset();
      for(int q=0; q < num_query; q++)
      {
         int offset = _data->subset[q].startIndex;
         int size = _data->subset[q].size;
         
         double ndcg = 0;
         eval_ndcg(offset, size, truncation, y_array, f_array, ndcg);
         mean_ndcg += ndcg;
         
         double auc = 0;
         eval_auc(offset, size, truncation, y_array, f_array, auc);
         mean_auc += auc;
         
      }
      mean_ndcg /= num_query;
      mean_ndcg *= 100;
      mean_auc /= num_query;
      mean_auc *= 100;
      
      printf("Average NDCG@%d: %.4f\n",truncation,mean_ndcg);
   }
}
