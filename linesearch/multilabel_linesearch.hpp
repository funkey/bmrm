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
 * Authors: Jin Yu (jin.yu@nicta.com.au)
 *          Choon Hui Teo (choonhui.teo@anu.edu.au)
 *
 * Created: (07/03/2008)
 *
 * Last Updated: 16/01/2009
 */

#ifndef _MULTILABEL_LINESEARCH_HPP_
#define _MULTILABEL_LINESEARCH_HPP_

#include "linesearch.hpp"
#include "multilabelvecdata.hpp"

#define MAX_INEXACT_LINESEARCH 20
#define MAX_LU_SEARCH 3

/**
 * Class for multilabel/multiclass line search
 */
class CMultilabel_Linesearch: public CLinesearch 
{
   protected:
      /**
       * Pointer to multilabelvecdata object
       */
      CMultilabelVecData* _data;

      /**
       * Vector of double-ended queue
       */
      std::vector<std::deque<lineWithLabel> > Ss;
      
      /**
       *  Vector of linear lines
       */
      std::vector<std::vector<Line> > Ls;
      
      /**
       * Pointer to a TheMatrix(1, numOfClass) object
       */
      TheMatrix* f_i;
      /**
       * Pointer to a TheMatrix(1, numOfClass) object
       */
      TheMatrix* df_i;
      
      /**
       * vector of iterators of the vector of active lines
       */
      std::vector<std::deque<lineWithLabel>::iterator> vec_iterator;
      
      /**
       * information about active linear lines such as active
       * labels, slope and bias etc.
       */
      std::vector<lineWithLabel> activeLine;
      
      
      /**
       * an object of a line
       */
      Line hLine;
      
      /**
       * Find all subdifferentiable/hinge points in [L, U],
       * and return the number of hinges
       *
       * @param f [read] prediction, i.e., XW
       * @param df [read] auxiliary vector XP
       * @param L [read] lower bound on the value the hinge
       * @param U [read] upper bound on the value the hinge
       * @param store_activeLabel [read] indicate whether to store active labels
       *
       */
      virtual int FindHinges(TheMatrix &f, TheMatrix &df, double L = 0.0,
                             double U = SML::INFTY, bool store_activelabel = true);
      
      /**
       * Given all linear lines (stored in Ls),
       * find all subdifferentiable/hinge points in [L, U],
       * and return the number of hinges.
       *
       * @param L [r/w] lower bound on the value the hinge
       * @param U [r/w] upper bound on the value the hinge
       * @param calculate_f [read] determine whether to evaluate
       * each linear function at L. If the function values at L have
       * been calculated and saved in line.value (line is saved in Ls),
       * set calculate_f = false.
       *
       *
       */
      virtual int FindHinges(double &L, double &U, bool calculate_f = true);
      
      /**
       * Store all linear lines
       *
       * @param f [read] prediction, i.e., XW
       * @param df [read] auxiliary vector XP
       *
       */
      virtual void ConstructLines(TheMatrix &f, TheMatrix &df);
      
      /**
       * Find [L, U] that contains the optimal step size,
       * by using gradient information.
       *
       * @param h [read] the value of lambda|p|^2
       * @param lambda_wp [read] lambda <w, p>
       * @param L [write] the lower bound
       * @param U [write] the upper bound
       *
       */
      virtual void FindLAndU(const double &h, const double &lambda_wp, double &L,
                             double &U);

      /**
       * Return the optimal step size.
       *
       * @param h [read] the value of lambda|p|^2
       * @param lambda_wp [read] lambda <w, p>
       * @param totoalLines [read] total number of linear lines in consideration
       * @param activeLine [write] active linear lines
       */
      virtual double FindOptStep(const double &h, const double &lambda_wp,
                                 int &totalLines, std::vector<lineWithLabel> &activeLine);
      
      /**
       * Return the optimal step size eta.
       *
       * @param h [read] the value of lambda|p|^2
       * @param lambda_wp [read] lambda <w, p>
       * @param totoalLines [read] total number of linear lines in consideration
       */
      virtual double FindOptStep(const double &h, const double &lambda_wp,
                                 int &totalLines);
      
      /**
       * Compute loss and supremum of subdifferential of sum of the pointwise maximum of
       * linear functions at a given step.
       * @param stepsize [read] the step size
       * @param sup_subgrad [write] the supremum of subdifferential
       * @param loss [write] loss
       */
      virtual void ComputeLossAndSupGrad(double &stepsize, double &sup_subgrad,
                                         double &loss);
      
      /**
       * Compute loss supremum of subdifferential of sum of the pointwise maximum of
       * linear functions at a given step at a given step in O(n log |Z|) using
       * hinge information stored in lists of sorted hinge stacks (Ss).
       * @param stepsize [read] the step size
       * @param sup_subgrad [write] the supremum of subdifferential
       * @param loss [write] loss
       */
      virtual void ComputeLossAndSupGradUsingSs(double &stepsize,
                                                double &sup_subgrad, double &loss);
      
      /**
       * Compute loss and supremum of subdifferential of sum of the pointwise maximum of
       * linear functions at a given step.
       * @param stepsize [read] the step size
       * @param sup_subgrad [write] the supremum of subdifferential
       * @param loss [write] loss
       * @param f [read] prediction, i.e., XW
       * @param df [read] auxiliary vector XP
       */
      virtual void ComputeLossAndSupGrad(double &stepsize, double &sup_subgrad,
                                         double &loss, TheMatrix &f, TheMatrix &df);
      
   public:
      
      /**
       * Class constructor.
       *
       * @param model [r/w] pointer to the model object
       * @param loss [r/w] pointer to the loss object
       * @param data [read] pointer to the data object
       */
      CMultilabel_Linesearch(CModel* model, CLoss* loss, CMultilabelVecData* data);
      
      /**
       * Class destructor
       */
      virtual ~CMultilabel_Linesearch() 
      {
         if (f_i) delete f_i;
         if (df_i) delete df_i;
      }
      
      /**
       * Execute the line search routine, and
       * Return the optimal step size.
       *
       * @param p [read] search direction
       * @param f [read] prediction, i.e., XW
       * @param df [read] auxiliary vector XP
       * @param lambda [read] regularizer constant
       * @param activeLine [write] active linear lines
       * @param calculate_loss [read] indicate whether to calculate loss
       */
      virtual double Linesearch(const TheMatrix &p, TheMatrix &f, TheMatrix &df,
                                const double &lambda, std::vector<lineWithLabel> &activeLine,
                                bool calculate_loss = false);
      
      /**
       * Execute the line search routine, and
       * Return the optimal step size.
       *
       * @param w [read] current iterate
       * @param p [read] search direction
       * @param f [read] prediction, i.e., XW
       * @param df [read] auxiliary vector XP
       * @param lambda [read] regularizer constant
       * @param calculate_loss [read] indicate whether to calculate loss
       * @param check_descent [read] indicate whether to check the given
       * direction p is a descent direction, e.g. the direction given by
       * BMRM bundle method solver may not be a descent direction at all, hence
       * there is no need to waste time in searching a step.
       */
      virtual double Linesearch(const TheMatrix &w, const TheMatrix &p,
                                TheMatrix &f, TheMatrix &df, const double &lambda,
                                bool calculate_loss = false, bool check_descent = false);
      
      
      /**
       * Execute the line search routine in two steps:
       * 1) find an interval [L, U] that contains the optimal
       * step size;
       * 2) run exact line search, and return optimal step size.
       *
       * @param p [read] search directioin
       * @param f [read] prediction, i.e., XW
       * @param df [read] auxiliary vector XP
       * @param lambda [read] regularizer constant
       * @param activeLine [write] active linear lines
       * @param calculate_loss [read] indicate whether to calculate loss
       */
      virtual double LinesearchInTwoSteps(const TheMatrix &p, TheMatrix &f,
                                          TheMatrix &df, const double &lambda,
                                          std::vector<lineWithLabel> &activeLine, bool calculate_loss = false);
      
      /**
       * Efficient bisection line search that returns a step size obeying
       * the sufficient decrease condition of the Wolfe conditions. This
       * line search calls FindHinge to identify active linear lines, and then caches
       * the slopes and offsets of those lines for efficient function value
       * evaluation.
       *
       * @param w [read] current iterate
       * @param p [read] search directioin
       * @param f [read] prediction, i.e., XW
       * @param df [read] auxiliary vector XP
       * @param lambda [read] regularizer constant
       * @param eta_0[read] Initial step size
       * @param c[read] Wolfe line search parameter
       * for sufficient decrease condition
       * @para a [read] step size shrinking factor
       */
      virtual double InexactLinesearch(const TheMatrix &w, const TheMatrix &p,
                                       TheMatrix &f, TheMatrix &df, const double &lambda, double eta_0 =
                                       1.0, double c = 1e-4, double a = 0.5);
      
};

#endif
