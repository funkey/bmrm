/* Copyright (c) 2009 NICTA
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
 * Created: (14/11/2007) 
 *
 * Last Updated: (21/12/2007)   
 */

#ifndef _MODEL_HPP_
#define _MODEL_HPP_

#include "sml.hpp"
#include <string>


/**  
 * A model is a container class which holds a weight vector.
 * A model is created in order to solve a regularized risk minimization problem. 
 * If a solution already exists (i.e. a model has already been created) 
 * then you can simply load the model and predict on test data. 
 *
 * If you need special methods for handling your data (e.g. reading in a
 * special format or writing out to a special format) then you need to
 * override the load and save methods. 
 * 
 */
class CModel 
{
   protected:
      /** Number of hyperplane
       */
      unsigned int _numOfW;
      
      /** Dimensionality of each hyperplane
       */
      unsigned int _dimOfW;
      
      /** Does the feature vector comes with additional feature for bias? 
       */
      bool _biasFeature;
      
      /** Weight vector
       */
      TheMatrix* _w;
      
      /** Initialize weight vector to zeros
       */ 
      TheMatrix* ZeroInit(const unsigned int& numOfW, const unsigned int& dimOfW, int type);
      
      /** Initialize weight vector to random uniform unit vector
       */
      TheMatrix* RandUniformInit(const unsigned int& numOfW, const unsigned int& dimOfW, int type);
      
      /** Initialize weight vector random vector from normal
       */
      TheMatrix* RandNormalInit(const unsigned int& numOfW, const unsigned int& dimOfW, int type);
      
   public:
      CModel() : _numOfW(0), _dimOfW(0), _biasFeature(0), _w(0) {}
      virtual ~CModel() { if(_w) delete _w; } 
      
      virtual void Reset(int type);
      /** Initialize the model
       *
       *  @param numOfW [read] number of hyperplane
       *  @param dimOfW [read] dimensionality of each hyperplane
       *  @param biasFeature [read] whether example comes with artificial bias feature
       */
      virtual void Initialize(const unsigned int &numOfW, const unsigned int &dimOfW, const bool &biasFeature);
      
      
      /** Initialize the model with values stored in the file
       *
       *  @param modelFilename [read] the file containing model information
       */
      virtual void Initialize(const std::string &modelFilename);
      
      
      /** Initialize the model with values stored in the file
       *
       *  @param modelFilename [read] the file containing model information
       *  @param dimOfW [read] the final dimensionality of w
       */
      virtual void Initialize(const std::string &modelFilename, const unsigned int& dimOfW);
      
      
      /** Whether model has been initialized
       */
      bool IsInitialized() { return (_dimOfW > 0 && _numOfW > 0); }
      
      
      /** Load the model from a file
       *
       *  @param modelFilename [read] Filename of the model
       */
      virtual void Load(const std::string &modelFilename);
      
      
      /** Save the model into a file
       *
       *  @param modelFilename [read] Filename of the model. If the model
       *     file name is set in the config, then fname is ignored.
       */
      virtual void Save(const std::string &modelFilename);
      
      
      /** Return a reference of the hyperplane(s)
       */
      TheMatrix& GetW() {return *_w;}
      
      
      /** Return the number of hyperplane
       */
      unsigned int GetNumOfW() {return _numOfW;}
      
      
      /** Return the dimensionality of hyperplane
       */
      unsigned int GetDimOfW() {return _dimOfW;}
};

#endif
