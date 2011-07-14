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
 * Created: (28/10/2008) 
 *
 * Last Updated:
 *
 */

#include <vector>
#include <fstream>
#include "configuration.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "solver.hpp"
#include "solverfactory.hpp"
#include "bmrm.hpp"
#include "loss.hpp"
#include "smmmulticlassloss.hpp"
#include "seqmulticlassdata.hpp"


using namespace std;

int main(int argc, char** argv)
{	
   // sanity check
   if(argc < 2) 
   {
      cout << "Usage: ./smmmc-bmrm-train config.file" << endl;
      cout << "Check the configfiles directory for examples" << endl;
      cout << "ERROR: No configuration file given!" << endl;
      exit(EXIT_FAILURE);
   }
    
 	               
   try {
      Configuration &config = Configuration::GetInstance();
      config.ReadFromFile(argv[1]);
      
      // DO NOT output decision function values and predicted labels to files after evaluation/prediction
      config.SetBool("Prediction.outputFvalAndLabels", false);
      
      CModel *model = new CModel();  	
      if(config.IsSet("Model.hotStartModel")) 
         model->Initialize(config.GetString("Model.hotStartModel"));		
      CSeqMulticlassData *data = new CSeqMulticlassData();		
      CLoss *loss = new CSMMMulticlassLoss(model, data);	
      //CBMRM* solver = new CBMRM;
      CSolver* solver = CSolverFactory::GetSolver(model,loss);

      solver->Train();		
      loss->Evaluate(model);
	
      string modelFn = "model";
      if(config.IsSet("Model.modelFile"))
         modelFn = config.GetString("Model.modelFile");
      model->Save(modelFn);
	
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
