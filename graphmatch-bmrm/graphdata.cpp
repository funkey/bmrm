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
 * Last Updated: 17/01/2009
 */

#ifndef _GRAPHDATA_CPP_
#define _GRAPHDATA_CPP_


#include "common.hpp"
#include "sml.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "timer.hpp"
#include "graphdata.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

//int WEIGHTS = 0;
//int NODES = 0;


adjmatrix::adjmatrix(int nodes)
    : NODES(nodes)
{
  mat = new int* [NODES];
  for (int i = 0; i < NODES; i ++)
  {
    mat[i] = new int [NODES];
    for (int j = 0; j < NODES; j ++) mat[i][j] = 0;
  }
}

adjmatrix::~adjmatrix()
{
  for (int i = 0; i < NODES; i ++)
    delete [] mat[i];
  delete [] mat;
}


CGraphData::CGraphData()
{
  Configuration &config = Configuration::GetInstance();
  string name;

  quad = config.GetInt("quadratic");
  int maxdata = config.GetInt("maxdata");

  NODES = 0;
  WEIGHTS = 0;

  NODES = config.GetInt("nnodes");
  WEIGHTS = config.GetInt("nfeatures");

  if(verbosity > 0)
     cout << "NODES=" << NODES << "  WEIGHTS=" << WEIGHTS << endl;

  assert(NODES > 0);
  assert(WEIGHTS > 0);
  assert (quad == 0 || quad == 1);

  impath = config.GetString("Data.path"); // path to examples
  name = config.GetString("Data.list");   // list of examples

  Ymem = new int [NODES];
  for (int i = 0; i < NODES; i ++)
    Ymem[i] = i;

  FILE* infile;
  infile = fopen(name.c_str(), "r");

  if (!infile)
  {
    cout << "Could not find " << name.c_str() << endl;
    exit(0);
  }

  char* imname = new char [256];
  char* cornername = new char [256];
  char* adjname = new char [256];
  char* scfname = new char [256];

  while (fscanf(infile, "%s %s %s %s", cornername, scfname, imname, adjname) == 4)
  {
     FILE* scffile = fopen((impath+"/"+string(scfname)).c_str(), "r");
    int** Gmem1 = new int* [NODES];

    //char* mname = new char [100];
    //memcpy(mname, imname, 100*sizeof(char));
    //imnames.push_back(mname);
    imnames.push_back(impath + "/" + string(imname));

    FILE* cornfile = fopen((impath+"/"+string(cornername)).c_str(), "r");
    if( cornfile == NULL)
    {
       ostringstream oss;
       oss << "Error: Could not open <" << string(cornername) << "> ";
       throw CBMRMException(oss.str(),"CGraphData()");
    }
    pair<double,double>* corns = new pair<double,double> [NODES];
    for (int i = 0; i < NODES; i ++)
    {
      float c1;
      float c2;
      fscanf(cornfile, "%f %f", &c1, &c2);
      corns[i].first = c1;
      corns[i].second = c2;
    }
    fclose(cornfile);
    corners.push_back(corns);

    for (int i = 0; i < NODES; i ++)
    {
      Gmem1[i] = new int [WEIGHTS];
      for (int j = 0; j < WEIGHTS; j ++)
      {
        int feat = 0;
        fscanf(scffile, "%d", &feat);
        Gmem1[i][j] = feat;
      }
    }
    Gmem.push_back(Gmem1);
    fclose(scffile);

    if (quad)
    {
      FILE* adjfile = fopen(adjname, "r");

      adjmatrix* Gadjmem1 = new adjmatrix(NODES);
      for (int i = 0; i < NODES; i ++)
      {
        for (int j = 0; j < NODES; j ++)
        {
          int adj = 0;
          fscanf(adjfile, "%d", &adj);
          Gadjmem1->mat[i][j] = adj;
        }
      }
      Gadjmem.push_back(Gadjmem1);
    }
  }

  delete [] imname;
  delete [] cornername;

  delete [] adjname;
  delete [] scfname;

  assert(quad == 0 || Gmem.size() == Gadjmem.size());

  if (config.GetString("Program.pairs") == "ALLPAIRS")
  {
    for (int i = 0; i < (int) Gmem.size(); i ++)
    {
      int j = 0;
      for (j = i+1; j < (int) Gmem.size(); j ++)
      {
        pairs.push_back(pair<int,int>(i,j));
        _G.push_back(Gmem[i]);
        _Gp.push_back(Gmem[j]);
        _Y.push_back(Ymem);
        if (quad)
        {
          _Gadj.push_back(Gadjmem[i]);
          _Gpadj.push_back(Gadjmem[j]);
        }
        if (maxdata > 0 && (int) _G.size() >= maxdata) break;
      }
      if (j < (int) Gmem.size()) break;
    }
  }
  else
  {
    for (int i = 0; i < (int) Gmem.size(); i += 2)
    {
      pairs.push_back(pair<int,int>(i,i+1));
      _G.push_back(Gmem[i]);
      _Gp.push_back(Gmem[i+1]);
      _Y.push_back(Ymem);
      if (quad)
      {
        _Gadj.push_back(Gadjmem[i]);
        _Gpadj.push_back(Gadjmem[i+1]);
      }
      if (maxdata > 0 && (int) _G.size() >= maxdata) break;
    }
  }

  _N = (int) _G.size();
}

CGraphData::~CGraphData()
{
  for (int i = 0; i < (int) Gmem.size(); i ++)
  {
    for (int j = 0; j < NODES; j ++)
      delete [] Gmem[i][j];
    delete [] Gmem[i];
    if (quad) delete Gadjmem[i];
  }

  //for (int i = 0; i < (int) imnames.size(); i ++) delete [] imnames[i];
  for (int i = 0; i < (int) corners.size(); i ++) delete [] corners[i];
}

#endif
