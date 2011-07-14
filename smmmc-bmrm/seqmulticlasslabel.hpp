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
 * Authors: Qinfeng Shi (qinfeng.shi@anu.edu.au)
 *          Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (14/04/2008) 
 *
 * Last Updated: (28/10/2008)
 */

#ifndef _SEQMULTICLASSLABEL_HPP_
#define _SEQMULTICLASSLABEL_HPP_

#include <vector>
#include <iostream>
#include <string>

#include "common.hpp"
#include "sml.hpp"

/** Class for semi markov model for sequence segmentation and classification
 *  
 *  Reference:
 *  1. Q. Shi , L. Wang, L. Cheng and A.J. Smola, 
 *     "Discriminative Human Action Segmentation and Recognition using Semi-Markov Model", 
 *     In IEEE CVPR 08, To Appear , 2008.
 *  2. Q. Shi, Y. Altun, A.J. Smola, and S.V.N. Vishwanathan,
 *     "Semi-Markov Models for Sequence Segmentation",
 *     in EMNLP, 2007.
 */


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
class CSeqMulticlassLabel
{        
public:
    CSeqMulticlassLabel();
    virtual ~CSeqMulticlassLabel() {}
    
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
    
    unsigned int getNumOfClass(){return numOfClass;}
    unsigned int getNumOfPerson(){return numOfPerson;}
    
    void Dump();
    
    
protected:
    /** Verbosity level
     */
    int seqmulticlasslabel_verbosity;
    
    /** labels
     */
    std::vector<seqlabel_struct> Y;
    
    /** Number of classes for segments
	One sequence has multiple class IDs on different segments.
    */
    unsigned int numOfClass;
    
    /** Number of classes for sequences
	In particular, it is called NumOfPerson to be consistent with ref[1].
	One sequence has a unique Person ID.
    */
    unsigned int numOfPerson;
    
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
