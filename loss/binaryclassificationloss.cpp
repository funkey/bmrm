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
 * Last Updated: (13/11/2007)   
 */


#include "binaryclassificationloss.hpp"

using namespace std;


/** Transform function value f := <w,x> into label
 *
 *  @param f [read/write] function values / labels
 */
void CBinaryClassificationLoss::Transform(TheMatrix& f)
{
	double* f_array = f.Data(); 
	int len = f.Length();
	for(int i=0; i < len; i++)
		f_array[i] = SML::sgn(f_array[i]);
}


/** Evaluate the performance of the model on _data
 *
 *  @param f [read] function values
 *  @param predict [write] prdicted labels
 */
void CBinaryClassificationLoss::Perf(TheMatrix& f, TheMatrix& predict)
{
 	double* p_array = predict.Data(); 
	double* y_array = _data->labels().Data();
	double* f_array = f.Data();
	unsigned int len = f.Length(); 
	
	// just checking
	assert(_data->labels().Length() == len);
	
	// evaluate performance using various metrics
	// compute TP, TN, FP, FN (i.e. contingency table)
	int TP = 0;    // true positive
	int TN = 0;    // true negative
	int FP = 0;    // false positive
	int FN = 0;    // false negative
	
	for(unsigned int i = 0; i < len; i++) 
	{
		if(p_array[i] > 0 && y_array[i] > 0)        TP++;
		else if(p_array[i] > 0 && y_array[i] <= 0)  FP++;
		else if(p_array[i] <= 0 && y_array[i] <= 0) TN++;
		else                                        FN++;
	}
	
	// performance measures
	double total = TP+TN+FP+FN;
	double precision = (TP>0)? 100.0*TP/(TP+FP) : 0.0;
	double recall = (TP>2) ? 100.0*TP/(TP+FN) : 0.0;
	double fone = ((precision + recall) > 0) ? 2.0*precision*recall/(precision+recall) : 0.0;
	
	// START : auc
	double auc = 0;
	vector<int> idx(len);
	for(size_t i=0; i < idx.size(); i++)
		idx[i] = i;
	
	indirect_less_than<double> ilt(f_array);
	sort(idx.begin(), idx.end(), ilt);  // ascending order
	
	// count swapped pairs i.e., misclassified points
	unsigned int n_pos = 0, n_neg = 0;
	double sum = 0;
	for(unsigned int i = 0; i < len; i++) 
	{
		if(y_array[idx[i]] < 0) 
		{
			sum += n_pos;
			n_neg++;
		}
		else
			n_pos++;
	}
	auc = 100.0 - 100.0*sum/n_neg/n_pos;
	// END
	
	// dump performanace to stdout
	cout << endl << "Performance:" << endl;
	cout << "1. Accuracy (%):   " << 100.0*(TP+TN)/total << endl;
	cout << "2. Precision (%):  " << precision << endl;
	cout << "3. Recall (%):     " << recall << endl;
	cout << "4. F1 score:       " << fone << endl;
	cout << "5. AUC:            " << auc << endl;
	cout << "6. True Positive:  " << TP << endl;
	cout << "7. True Negative:  " << TN << endl;
	cout << "8. False Positive: " << FP << endl;
	cout << "9. False Negative: " << FN << endl;
}
