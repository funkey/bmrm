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
 * Last Updated: (26/01/2008)   
 */

#ifndef _DATA_HPP_
#define _DATA_HPP_

#include <iostream>
#include <vector>


/** Base class for encapsulating a data container.
 *  This class acts almost like a void*
 */
class CData 
{
   protected:
      /** verbosity
       */
      int verbosity;     
              
   public:
      CData();      
      virtual ~CData() {};    
      
      /** 
       * Return the number of examples in this sub-dataset.
       */
      virtual unsigned int slice_size(void) const = 0;
      
      /** 
       * Return the total number of examples in this dataset.
       */
      virtual unsigned int size(void) const = 0;
      
      /** 
       * Return the dimension of example in this dataset.
       */
      virtual unsigned int dim(void) const = 0;
      
      /** 
       * Return true if we use a bias feature. False otherwise. 
       */
      virtual bool bias(void) const = 0;      
};

#endif
