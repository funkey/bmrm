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

#ifndef _NDCGRANKLOSS_HPP_
#define _NDCGRANKLOSS_HPP_

#include <vector>
#include "sml.hpp"
#include "rankloss.hpp"
#include "model.hpp"

/** NDCG Ranking loss.
 *  Will create a superclass such as CRankLoss for all ranking loss soon...
 */
class CNDCGRankLoss : public CRankLoss
{
   protected:               
      int truncation;
      std::vector<int> current_ideal_pi;
      std::vector< std::vector<int> > sort_vectors; /* sort indices per query */
      int max_subset_size;
      std::vector<double> c;
      int *pi;
      std::vector < std::vector<double> > bs;
      std::vector<double> a;
      
      void Loss(double& loss, TheMatrix& f);
      void LossAndGrad(double& loss, TheMatrix& f, TheMatrix& l);
      void delta(int size, std::vector<double> a, std::vector<double> b, 
                 int *pi, double &value);
      void find_permutation(int size, int offset, std::vector<double> a, 
                            std::vector<double> b, std::vector<double> c, 
                            double *f, int *pi);
      void compute_coefficients(int offset, int size, double *y, 
                                std::vector<int> ideal_pi, 
                                std::vector<double> a, std::vector<double> b);
      double get(double *f, int offset, int i);
      void add(TheMatrix &l, int offset, int i, double value);
      void compute_b(int offset, int size, double *y, 
                     std::vector<int> ideal_pi, 
                     std::vector<double> &input_b);
      void compute_a(int size, std::vector<double> &input_a);
      
   public:    
      
      CNDCGRankLoss(CModel* &model, CVecData* &data);
      virtual ~CNDCGRankLoss() {if (pi){free(pi);pi = 0;}}
};

#endif
