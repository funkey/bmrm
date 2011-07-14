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
 * Authors: Simon Guenter
 *          Jin Yu (jin.yu@nicta.com.au)
 *          Choon Hui Teo (choonhui.teo@anu.edu.au)
 *
 * Created: (14/02/2008)
 *
 * Last Updated: 16/01/2009
 */

#ifndef _HINGE_LINESEARCH_HPP_
#define _HINGE_LINESEARCH_HPP_

#include <algorithm>
#include "linesearch.hpp"

class CHinge_Linesearch: public CLinesearch 
{
   protected:
      std::vector<ScalarWithIndex> alpha;
      std::vector<size_t> hinges;
      
   public:
      /**
       * Class constructor.
       *
       * @param model [read] pointer to the model object
       * @param loss [read] pointer to the loss object
       * @param data [read] pointer to the data object
       */
      CHinge_Linesearch(CModel* model, CLoss* loss, CData* data);
      
      /**
       * Destructor
       */
      virtual ~CHinge_Linesearch() {}
      
      /**
       * Execute the line search routine, and
       * Return step size eta.
       *
       * @param p [read] search directioin
       * @param f [read] prediction, i.e., y.(Xw)
       * @param df [read] auxiliary vector y.(Xp)
       * @param lambda [read] regularizer constant
       * @param delta [write] indicator vector: 1:error, 0:correct
       */
      virtual double Linesearch(const TheMatrix &p, TheMatrix &f, TheMatrix &df,
                                const double &lambda, TheMatrix &delta);
      
      /**
       * Return a list hinge indices.
       */
      virtual std::vector<size_t> GetHinges(void) const 
      {
         return hinges;
      }
};

#endif
