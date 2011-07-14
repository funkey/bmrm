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
 * Created: 29/05/2008
 *
 * Last Updated:
 *
 */

#include <vector>
#include <fstream>
#include "configuration.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "modelfactory.hpp"
#include "bmrm.hpp"
#include "loss.hpp"
#include "linesearchlossfactory.hpp"
#include "data.hpp"
#include "datafactory.hpp"


using namespace std;

int main(int argc, char** argv)
{	
   // sanity check
   if(argc < 2) 
   {
      cout << "Usage: ./ls-bmrm-train config.file" << endl;
      cout << "Check the configfiles directory for examples" << endl;
      cout << "ERROR: No configuration file given!" << endl;
      exit(EXIT_FAILURE);
   }
	
   CData* data = 0;
   CLoss* loss = 0;
   CBMRM* solver = 0;
   CModel* model = 0;
   
   try {
      Configuration &config = Configuration::GetInstance();
      config.ReadFromFile(argv[1]);
      
      // DO NOT output decision function values and predicted labels to files after evaluation
      config.SetBool("Prediction.outputFvalAndLabels", false);

      model = CModelFactory::GetModel();  
		
      if(config.IsSet("Model.hotStartModel")) 
         model->Initialize(config.GetString("Model.hotStartModel"));		
      data = CDataFactory::GetData();		
      loss = CLinesearchLossFactory::GetLoss(model, data); // loss will initialize model if model is not hotstarted                                             
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
            solver = new CBMRM(model, loss);
            solver->Train();		
            loss->Evaluate(model);
			
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
      else
      {
         solver = new CBMRM(model, loss);
         solver->Train();		
         loss->Evaluate(model);
         
         string modelFn = "model";
         if(config.IsSet("Model.modelFile"))
            modelFn = config.GetString("Model.modelFile");
         model->Save(modelFn);
      }                
      
      // cleaning up 
      if(solver) delete solver;
      if(model) delete model;
      if(loss) delete loss;
      if(data) delete data;		
   } 
   catch(CBMRMException e) 
   {
      cerr << e.Report() << endl;
   }
   
   return EXIT_SUCCESS;	
}
