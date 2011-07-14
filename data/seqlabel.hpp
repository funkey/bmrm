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
 * Created: (28/01/2008) 
 *
 * Last Updated:
 */

#ifndef _SEQLABEL_HPP_
#define _SEQLABEL_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"


/** Container for sequence labels. Mainly for the application where
 *  sequence has a particular id, and the sequence element has an associated class.
 *  
 *   Label file format: rows of <label>
 *   where
 *   <label> .=. <pid>:<N> <pos>:<type> [<pos>:<type>...]
 *   <pid>   .=. sequence id (natural number)
 *   <N>     .=. natural number
 *   <pos>   .=. element of sequence (natural number)  
 *   <type>  .=. class/type associated to a segment of sequence (natural number)
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 *   2.  For centralized dataset and serial computation mode only!
 */
class CSeqLabel
{        
   public:
      CSeqLabel();
      virtual ~CSeqLabel() {}
      
      /** label type
         */
      struct seqlabel_struct {
            unsigned int ID;
            std::vector<unsigned int> pos;
            std::vector<unsigned int> type;
      };
      
      
      /** Return the labels. 
       */
      virtual const std::vector<seqlabel_struct>& labels(void){ return Y; }         
      
        void Dump();
      /**  subsets
       */
      //vector<Subset> subset;    
      
   protected:
      /** Verbosity level
       */
      int seqlabel_verbosity;
        
      
      /** labels
       */
      std::vector<seqlabel_struct> Y;
      
      /** Structure to keep some info. of subsets of examples in the dataset
       */
      //struct Subset {
      //        int ID;
      //        int startIndex;
      //        int size;
      //};                   
      
      /** Number of labels in the sub-dataset
       */
      unsigned int numOfLabel;
      
      /** Number of labels in the WHOLE dataset
       *  This is the same as numOfLabel in centralised dataset, serial computation mode
       */
      unsigned int numOfAllLabel;
      
      /** The first example of the dataset (skip #startExample# examples before reading the rest)
       */
      unsigned int startExample;
      
      /** Name of the file containing label vectors
       */
      std::string labelFile;
      
      /** Allocate label vector (or matrix) and load label from label file
       */
      virtual void LoadLabels(); 
      
};

#endif
