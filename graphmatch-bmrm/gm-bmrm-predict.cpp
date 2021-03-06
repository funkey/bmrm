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
 * Created: (27/10/2008) 
 *
 * Last Updated:
 */

#include "configuration.hpp"
#include "bmrm.hpp"
#include "sml.hpp"
#include "graphmatchloss.hpp"
#include "graphdata.hpp"
#include "model.hpp"

using namespace std;


int main(int argc, char** argv)
{    
   // sanity check
   if(argc < 2) 
   {
      cout << "Usage: gm-bmrm-predict config.file" << endl;
      cout << "Check the configfiles directory for examples" << endl;
      cout << "ERROR: No configuration file given!" << endl;
      exit(EXIT_FAILURE);
   }
   
   try {
      Configuration &config = Configuration::GetInstance();
      config.ReadFromFile(argv[1]);
      
      // output decision function values and predicted labels to files after evaluation/prediction
      config.SetBool("Prediction.outputFvalAndLabels", true);
      
      string modelFilename = config.GetString("Model.modelFile");
      
      CGraphData *data = new CGraphData();     
      CModel *model = new CModel();
      model->Initialize(modelFilename, data->dim());
      CGraphMatchLoss *loss = new CGraphMatchLoss(model,data);
      loss->Evaluate(model);
       
      if (model )delete model;
      if(loss) delete loss;
      if(data) delete data;
   }
   catch(CBMRMException e) {
      cout << e.Report() << endl;
   }
   
   return EXIT_SUCCESS;
}
