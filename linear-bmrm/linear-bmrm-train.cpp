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
 * Created: (02/11/2007) 
 *
 * Last Updated: 15/01/2009
 *
 * Note: [chteo:100108:1743] now with support for model selection over lambda parameters.
 *       The model files will be appended with "_lambda_XXX" where XXX is the lambda used 
 *       to produce that model.
 */

#include <vector>
#include <fstream>
#include "configuration.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "modelfactory.hpp"
#include "solver.hpp"
#include "solverfactory.hpp"
#include "bmrm.hpp"
#include "loss.hpp"
#include "lossfactory.hpp"
#include "data.hpp"
#include "datafactory.hpp"

#ifdef PARALLEL_BMRM
#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR
#include "mpi.h"
#endif

using namespace std;

int main(int argc, char** argv)
{	
   // sanity check
   if(argc < 2) 
   {
      cout << "Usage: ./linear-bmrm-train config.file" << endl;
      cout << "Check the configfiles directory for examples" << endl;
      cout << "ERROR: No configuration file given!" << endl;
      exit(EXIT_FAILURE);
   }
   
#ifdef PARALLEL_BMRM   
   int myProcID = -1;
   int nProc = -1;
   
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &myProcID);
   MPI_Comm_size(MPI_COMM_WORLD, &nProc);
   
   MASTER(myProcID)
   {
      cout << nProc << " machines found!" << endl;
   }
#endif

   CData* data = 0;
   CLoss* loss = 0;
   CSolver* solver = 0;
   CModel* model = 0;
   
   try {
      Configuration &config = Configuration::GetInstance();
      config.ReadFromFile(argv[1]);
      
      // DO NOT output decision function values and predicted labels to files after evaluation
      config.SetBool("Prediction.outputFvalAndLabels", false);
      
      model = CModelFactory::GetModel();  
      if(config.IsSet("Model.hotStartModel")) 
         model->Initialize(config.GetString("Model.hotStartModel"));		
#ifdef PARALLEL_BMRM
      data = CDataFactory::GetData(myProcID,nProc);
#else
      data = CDataFactory::GetData();
#endif

      loss = CLossFactory::GetLoss(model, data); // loss will initialize model if model is not hot-started

      vector<double> lambdas;
      if(config.IsSet("BMRM.lambdas"))
      {
         lambdas = config.GetDoubleVector("BMRM.lambdas");
         cout << "Main(): Multiple lambda parameters detected! Training with each of them..." << endl;
         
         for(size_t i=0; i < lambdas.size(); i++)
         {
            config.SetDouble("BMRM.lambda", lambdas[i]);
            cout << "\n[Learning using lambda: " << lambdas[i] << "]" << endl;
            
            if(solver) delete solver;
#ifdef PARALLEL_BMRM
            solver = CSolverFactory::GetSolver(model,loss,myProcID,nProc);
#else
            solver = CSolverFactory::GetSolver(model,loss);
#endif
            solver->Train();		
#ifndef PARALLEL_BMRM
            // in parallel computation, master holds only a subset of whole dataset
            // so evaluation of training will be bogus
            loss->Evaluate(model);
#endif
#ifdef PARALLEL_BMRM
            MASTER(myProcID)
#endif
            {
               ostringstream oss;
               oss << lambdas[i];
               string lambda_str = oss.str();
               string modelFn = "model_lambda_" + lambda_str;
               if(config.IsSet("Model.modelFile"))
                  modelFn = config.GetString("Model.modelFile");
               modelFn = modelFn + "_lambda_" + lambda_str;
               model->Save(modelFn);
            }
         }
      }
      else
      {
#ifdef PARALLEL_BMRM
            solver = CSolverFactory::GetSolver(model,loss,myProcID,nProc);
#else
            solver = CSolverFactory::GetSolver(model,loss);
#endif
         solver->Train();		
#ifdef PARALLEL_BMRM
         MASTER(myProcID)
#endif
         {
#ifndef PARALLEL_BMRM
            // in parallel computation, master holds only a subset of whole dataset
            // so evaluation of training will be bogus
            loss->Evaluate(model);
#endif
            string modelFn = "model";
            if(config.IsSet("Model.modelFile"))
               modelFn = config.GetString("Model.modelFile");
            model->Save(modelFn);
         }
      }                
      
      // cleaning up 
      if(solver) delete solver;
      if(model) delete model;
      if(loss) delete loss;
      if(data) delete data;		
   } 
   catch(CBMRMException e) {
      cerr << e.Report() << endl;
   }
   
   return EXIT_SUCCESS;	
}
