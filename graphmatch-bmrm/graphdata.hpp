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
 * Authors: Julian McAucley (julian.mcaculey@nicta.com.au)
 *
 * Created: (27/10/2008) 
 *
 * Last Updated:
 */

#ifndef _GRAPHDATA_HPP_
#define _GRAPHDATA_HPP_

#include "common.hpp"
#include "data.hpp"
#include <vector>
#include <iostream>

//extern int WEIGHTS;
//extern int NODES;

class adjmatrix
{
public:
   adjmatrix(int nodes);
    ~adjmatrix();
    
   int** mat;
   int NODES;
};

class CGraphData: public CData
{
   public:
      // constructor
      CGraphData();
      virtual ~CGraphData();
      
      std::vector<int**> _G;
      std::vector<int**> _Gp;
      std::vector<int*> _Y;
      
      std::vector<adjmatrix*> _Gadj;
      std::vector<adjmatrix*> _Gpadj;
      
      std::vector<std::pair<int,int> > pairs;
      
      int WEIGHTS;
      int NODES;
      
      int quad;
      int _N;
      
      unsigned int numOfExample;
      unsigned int numOfAllExample;
      
      std::string impath;
      //std::vector<char*> imnames;
      std::vector<std::string> imnames;
      std::vector<std::pair<double,double>*> corners;
      
      virtual bool bias(void) const { return false; }
      virtual unsigned int dim(void) const { return WEIGHTS + quad; }
      
      virtual unsigned int slice_size(void) const { return _N; }
      virtual unsigned int size(void) const { return _N; }
      
   protected:
      std::vector<int**> Gmem;
      int* Ymem;
      std::vector<adjmatrix*> Gadjmem;
      
   private:
      CGraphData(const CGraphData&);
};

#endif
