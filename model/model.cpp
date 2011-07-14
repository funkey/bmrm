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
 * Last Updated: (14/11/2007)   
 */

#include <fstream>
#include <sstream>

#include "common.hpp"
#include "model.hpp"
#include "configuration.hpp"

using namespace std;

TheMatrix* CModel::ZeroInit(const unsigned int& numOfW, const unsigned int& dimOfW, int type)
{
   TheMatrix *tmp = new TheMatrix(1, numOfW*dimOfW, type);
   return tmp;
}


void CModel::Reset(int type)
{
   _w = ZeroInit(1, _numOfW*_dimOfW, type);
}


TheMatrix* CModel::RandUniformInit(const unsigned int& numOfW, const unsigned int& dimOfW, int type)
{
   TheMatrix *tmp = new TheMatrix(1, numOfW*dimOfW, type);
   srand(time(0));
   long idum = -rand();
   for(unsigned int i=0; i< numOfW*dimOfW; i++)
      tmp->Set(i, ran1(&idum));
   return tmp;
}


TheMatrix* CModel::RandNormalInit(const unsigned int& numOfW, const unsigned int& dimOfW, int type)
{
   TheMatrix *tmp = new TheMatrix(1, numOfW*dimOfW, type);
   srand(time(0));
   long idum = -rand();
   for(unsigned int i=0; i< numOfW*dimOfW; i++)
      tmp->Set(i, gasdev(&idum));
   return tmp;
}


void CModel::Initialize(const unsigned int &numOfW, const unsigned int &dimOfW, const bool &biasFeature)
{
   assert(! IsInitialized());
   
   _numOfW = numOfW;
   _dimOfW = dimOfW;
   _biasFeature = biasFeature;
   
   Configuration &config = Configuration::GetInstance();
   string wInitType = "ZERO";
   if(config.IsSet("Model.initialization"))
   {
      wInitType = config.GetString("Model.initialization");
      if(wInitType == "RANDOM_UNIFORM")
      {
         _w = RandUniformInit(1, _numOfW*_dimOfW, SML::DENSE);
         printf("CModel:: random uniformly distributed weight vector initialized\n");
      }        
      else if(wInitType == "RANDOM_NORMAL")
      {
         _w = RandNormalInit(1, _numOfW*_dimOfW, SML::DENSE);
         printf("CModel:: random normaly distributed weight vector initialized\n");
      }
   }
   else
      _w = ZeroInit(1, _numOfW*_dimOfW, SML::DENSE);
}


void CModel::Initialize(const string &fname)
{
   ifstream fp(fname.c_str());
   if(! fp.good()) 
   {
      throw CBMRMException("Cannot open file <" + fname + ">\n","CModel::Initialize()");
   }   
   
   string line, dummy;
   istringstream iss;
   
   getline(fp,line);
   iss.str(line);
   iss >> dummy >> _numOfW >> _dimOfW;
   
   getline(fp,line);
   iss.str(line);
   iss >> dummy >> _biasFeature;
   
   _w = new TheMatrix(1, _numOfW * _dimOfW, SML::DENSE);
   
   double value = 0.0;
   for(unsigned int i=0; i<_numOfW * _dimOfW; i++) 
   {
      fp >> value;
      _w->Set(i,value);
   }
   
   fp.close();
}


void CModel::Initialize(const string &fname, const unsigned int &dimOfW)
{
   ifstream fp(fname.c_str());
   if(! fp.good()) 
   {
      throw CBMRMException("Cannot open file <" + fname + ">\n","CModel::Initialize()");
   }   
   
   string line, dummy1, dummy2;
   istringstream iss;
   
   _dimOfW = dimOfW;
   
   // get numOfW
   getline(fp,line);
   iss.str(line);
   iss >> dummy1 >> _numOfW >> dummy2;
   
   // get biasFeature flag
   getline(fp,line);
   iss.str(line);
   iss >> dummy1 >> _biasFeature;
   
   _w = new TheMatrix(1, _numOfW * _dimOfW, SML::DENSE);
   
   fp.close();
   
   Load(fname);
}


void CModel::Load(const string &fname)
{
   if(_w == 0)
   {
      CBMRMException("weight vector has not been initialized yet!","CModel::Load()");
   }      
   ifstream fp(fname.c_str());
   if(! fp.good()) 
   {
      throw CBMRMException("Cannot open file <" + fname + ">\n","CModel::Load()");
      return;
   }   
   
   unsigned int in_numOfW = 0, in_dimOfW = 0, in_bias = 0;
   string line, dummy;
   istringstream iss;
   
   getline(fp,line);
   iss.str(line);
   iss >> dummy >> in_numOfW >> in_dimOfW;
   
   getline(fp,line);
   iss.str(line);
   iss >> dummy >> in_bias;
   
   int dimFeature = (_biasFeature) ? _dimOfW-1 : _dimOfW;
   int dimWeight = in_dimOfW - in_bias;
   int truncation = max(0, dimWeight - dimFeature);
   int extension = max(0, dimFeature - dimWeight);
   
   
   if(in_numOfW != _numOfW) 
   {
      ostringstream oss;
      oss << "number of W in model (" << in_numOfW 
          << ")  does not match with the number of W in data (" 
          << _numOfW << ")";
      throw CBMRMException(oss.str(),"CModel::Load()");
   }
   
   if(truncation > 0) 
   {
      cout << "Warning: dim(feature vector)=" << dimFeature 
           << " > dim(weight vector)=" << dimWeight 
           << ". Last " << truncation 
           << " entries of weight vector will be truncated!" << endl;
   }
   else if(extension > 0) 
   {
      cout << "Warning: dim(feature vector)=" << dimFeature 
           << " > dim(weight vector)=" << dimWeight 
           << ". Weight vector will be extended (with " 
           << extension << " zeros)!" << endl;
   }
   
   double value = 0.0;
   for(unsigned int i=0; i<_numOfW; i++) 
   {
      unsigned int j = 0;
      // load weight vector up to truncation            
      for(j=0; j < in_dimOfW - in_bias - truncation; j++)             
      {
         fp >> value;
         _w->Set(i*_dimOfW + j,value);
      }
      
      // get rid of the truncated entries from file
      for(int k=0; k<truncation; k++)
         fp >> value;
      
      // extension of weight vector
      for(int ext_i=0; ext_i<extension; ext_i++, j++) 
         _w->Set(i*_dimOfW + j, 0.0);
      
      // bias feature of weight vector
      if((in_bias == 1) && _biasFeature) 
      {
         // j has been increased in previous for loop
         fp >> value;
         _w->Set(i*_dimOfW + j,value);         
      }
   }
   
   if((in_bias == 1) && (not _biasFeature)) 
      cout << "Warning: Feature vectors do not come with bias features. "
           << "The bias feature in weight vector will be dropped!" 
           << endl;   
   
   fp.close();
}


void CModel::Save(const string &fname)
{
   ofstream outfp(fname.c_str());
   double value = 0.0;
   
   if(! outfp.good())
   {
      ostringstream oss;
      oss << "failed to open file <" << fname << "> for saving model";
      throw CBMRMException(oss.str(),"CModel::Save()");      
   }
   
   
   // dimensionality of W and information about bias feature
   outfp << "# " << _numOfW << " " << _dimOfW 
         << " (number of hyperplane, dimensionality of hyperplane)" << endl;
   outfp << "# " << _biasFeature << " (added bias feature? 0:no 1:yes)" << endl;
   
   for(unsigned int i=0; i < _numOfW * _dimOfW; i++) 
   {
      _w->Get(i, value);
      outfp << value << endl;
   }
   
   outfp.close();  
}
