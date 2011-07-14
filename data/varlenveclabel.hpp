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
 * Created: (26/01/2008) 
 *
 * Last Updated:
 */

#ifndef _VARLENVECLABEL_HPP_
#define _VARLENVECLABEL_HPP_

#include <vector>
#include <iostream>

#include "common.hpp"
#include "sml.hpp"


/** Container for variable length labels 
 *    
 *   Label file format: rows of <label>
 *   where
 *   <label>   .=. <value_1>,[value_2>,...,<value_k>]
 *   <value_i> .=. i-th label (scalar value) of a particular example
 *
 * Notes:
 *   1.  The data loading starts during the object construction.
 */
class CVarLenVecLabel
{       
   protected:
      /** Verbosity level
       */
      int varlenveclabel_verbosity;
      
      /** Label vectors
       */
      std::vector<std::vector<double> > Y;
      
      /** Average length of labels
       */
      double avgLabelDimension;
      
      /** Smallest label (needed in determining the number of classes in multi-class classification)
       */
      double minLabel;
      
      /** Largest label (needed in determining the number of classes in multi-class classification)
       */
      double maxLabel;
      
      /** Number of labels in the sub-dataset
       */
        unsigned int numOfLabel;
      
      /** Number of labels in the WHOLE dataset
       *  This is the same as numOfLabel in centralised dataset, serial computation mode
       */
      unsigned int numOfAllLabel;
      
      /** The index of the first label of the dataset (skip #startLabelIndex# examples before reading the rest)
       */
      unsigned int startLabelIndex;
      
      /** Name of the file containing label vectors
       */
      std::string labelFile;
      
      /** A string template for floating point value
       */
      std::string scalar_value_format;

      /** User-provided valued to be added to each label
       */
      double labelIncrement;

      
      /** Scan label file to determine some label set properties such as type of label
       */
      virtual void ScanLabelFile();
      
      /** Allocate label vector (or matrix) and loat label from label file
       */
      virtual void LoadLabels(); 
      
   public:
      CVarLenVecLabel();
      virtual ~CVarLenVecLabel();
      
      /** Load labels.
       *  Load only "nlabels" labels starting from "startlabel"
       */
      virtual void LoadVarLenLabel(unsigned int startlabel, unsigned int nlabels);
        
      /** Return the labels. 
       */
      virtual const std::vector<std::vector<double> >& labels(void){ return Y; }
      
      
      /** 
       * Return the maximum of labels
       */
      double MaxLabel() const {return maxLabel;}
      
      
      /** 
       * Return the minimum of labels
       */
      double MinLabel() const {return minLabel;}
      
      /** Show label
       */
      virtual void ShowLabel();
};

#endif
