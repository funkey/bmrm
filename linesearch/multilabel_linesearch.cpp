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
 * Authors: Jin Yu (jin.yu@nicta.com.au)
 *          Choon Hui Teo (choonhui.teo@anu.edu.au)
 *
 * Created: (07/03/2008)
 *
 * Last Updated: 16/01/2009
 */

#include "multilabel_linesearch.hpp"

using namespace std;

/**
 * A comparison operation that is called by  the 'sort' operation
 */
inline bool descent(const Line &i, const Line &j) 
{
   return (i.value > j.value);
}

/**
 * A comparison operation that is called by  the 'lower_bound' operation
 */
inline bool comp(const lineWithLabel &i, const lineWithLabel &j) 
{
   return (i.step > j.step);
}

CMultilabel_Linesearch::CMultilabel_Linesearch(CModel* model, CLoss* loss, CMultilabelVecData* data) 
   : CLinesearch(model, loss, data), 
     _data(data),
     f_i(0), 
     df_i(0) 
{
   f_i = new TheMatrix();
   df_i = new TheMatrix();
   vec_iterator.clear();
   Ss.clear();
   Ls.clear();
   activeLine.clear();
   hLine = Line(0.0, LabelPair(-1, -1), 0.0, 0.0);
   
}

void CMultilabel_Linesearch::ComputeLossAndSupGrad(double &stepsize, double &sup_subgrad, double &loss) 
{
   sup_subgrad = 0.0;
   loss = 0.0;
   for (vector<vector<Line> >::iterator it = Ls.begin(); it < Ls.end(); it++) {
      double max_slope = -SML::INFTY, max_value = -SML::INFTY;
      for (vector<Line>::iterator itt = it->begin(); itt < it->end(); itt++) {
         itt->value = itt->bias + stepsize * itt->slope;
         
         if (itt->value > max_value || (itt->value == max_value
                                        && itt->slope > max_slope)) {
            max_slope = itt->slope;
            max_value = itt->value;
            
         }
      }
      
      sup_subgrad += max_slope;
      loss += max(0.0, max_value);
   }
}

void CMultilabel_Linesearch::ComputeLossAndSupGradUsingSs(double &stepsize,
                                                          double &sup_subgrad, double &loss) {
   loss = 0.0;
   sup_subgrad = 0.0;
   lineWithLabel line;
   line.step = stepsize;
   
   for (vector<deque<lineWithLabel> >::iterator it = Ss.begin(); it < Ss.end(); it++) {
      deque<lineWithLabel>::iterator up = upper_bound(it->begin(), it->end(),
                                                      line, comp);
      if (up->step >= stepsize) {
         up--;
      }
      double max_value = up->bias + stepsize * up->slope;
      loss += max(0.0, max_value);
      sup_subgrad += up->slope;
   }
}

void CMultilabel_Linesearch::ComputeLossAndSupGrad(double &stepsize,
                                                   double &sup_subgrad, 
                                                   double &loss, 
                                                   TheMatrix &f, 
                                                   TheMatrix &df) 
{
   
   const vector<vector<double> >& Y = _data->labels();
   double slope = 0.0, bias = 0.0, value = 0.0;
   loss = 0.0;
   sup_subgrad = 0.0;
   for (int i = 0; i < numElements; i++) {
      f.CreateMatrixRowView(i, f_i);
      df.CreateMatrixRowView(i, df_i);
      double max_slope = 0.0, max_value = 0.0;
      
      //create linear lines, and set 'Line.value' to be the value
      //evaluated at 0
      for (int k = 0; k < numClass; k++) {
         bool isTrueLabel = false;
         //check whether k is in the correct label set
         for (uint j = 0; j < Y[i].size(); j++) {
            if (k + 1 == Y[i][j]) {
               isTrueLabel = true;
               break;
            }
         }
         
         //if k is not in the correct label set
         //then add numCorrectLabel number of linear lines:
         //l := bias + stepsize*slope, where the values of bias and slope
         //depend on the pair of label (k, z), where z
         //is in the correct label set
         if (!isTrueLabel) {
            for (uint j = 0; j < Y[i].size(); j++) {
               double f_k, df_k, f_y, df_y;
               f_i->Get(k, f_k);
               df_i->Get(k, df_k);
               f_i->Get(int(Y[i][j]) - 1, f_y);
               df_i->Get(int(Y[i][j]) - 1, df_y);
               //bias
               bias = (1.0 + f_k - f_y);
               //slope
               slope = df_k - df_y;
               //value at lower bound
               value = bias + stepsize * slope;
               
               if (value > max_value || (value == max_value && slope
                                         > max_slope)) {
                  max_slope = slope;
                  max_value = value;
                  
               }
               
            }
         }
      }
      
      sup_subgrad += max_slope;
      loss += max(0.0, max_value);
   }
}

void CMultilabel_Linesearch::ConstructLines(TheMatrix &f, TheMatrix &df) {
   
   const vector<vector<double> >& Y = _data->labels();
   Line oneLine;
   vector<Line> lines;
   vector<Line>::iterator it_lines;
   //a horizontal line, y = 0. note we set trueLabel=-1
   //and label=-1 since this line is not associated with any label.
   //Line hLine = Line(0.0, LabelPair(-1, -1), 0.0, 0.0);
   for (int i = 0; i < numElements; i++) {
      
      f.CreateMatrixRowView(i, f_i);
      df.CreateMatrixRowView(i, df_i);
      lines.clear();
      if (Ls.size() < uint(numElements)) {
         //add a horizontal line for max(0,.) operation of hinge loss
         lines.push_back(hLine);
      } else {
         it_lines = (Ls.at(i)).begin();
         *(it_lines) = hLine;
         it_lines++;
      }
      //create linear lines, and set 'Line.value' to be the value
      //of the linear line evaluated at 0
      for (int k = 0; k < numClass; k++) {
         bool isTrueLabel = false;
         //check whether k is in the correct label set
         for (uint j = 0; j < Y[i].size(); j++) {
            if (k + 1 == Y[i][j]) {
               isTrueLabel = true;
               break;
            }
         }

         //if k is not in the correct label set,
         //then add numCorrectLabel number of linear lines:
         //l := bias + stepsize*slope, where the values of bias and slope
         //depend on the pair of label (k, z), where z
         //is in the correct label set
         if (!isTrueLabel) {
            for (uint j = 0; j < Y[i].size(); j++) {
               double f_k, df_k, f_y, df_y;
               f_i->Get(k, f_k);
               df_i->Get(k, df_k);
               f_i->Get(int(Y[i][j]) - 1, f_y);
               df_i->Get(int(Y[i][j]) - 1, df_y);
               //bias
               oneLine.bias = (1.0 + f_k - f_y);
               //slope
               oneLine.slope = df_k - df_y;
               //associated label pair. note k is shifted by 1
               oneLine.labels = LabelPair(int(Y[i][j]), k + 1);
               
               //this is done to prevent numerical precision error
               oneLine.value *= (abs(oneLine.value) > MYZERO);
               oneLine.slope *= (abs(oneLine.slope) > MYZERO);
               oneLine.bias *= (abs(oneLine.bias) > MYZERO);
               if (Ls.size() < uint(numElements)) {
                  lines.push_back(oneLine);
               } else {
                  *(it_lines) = oneLine;
                  it_lines++;
               }
            }
         }
      }

      if (Ls.size() < uint(numElements)) {
         Ls.push_back(lines);
      }
      
   } // end of creating all the linear lines
   
}

void CMultilabel_Linesearch::FindLAndU(const double &h,
                                       const double &lambda_wp, double &L, double &U) {
   
   double sup_subgrad = 0.0;
   double test_point = 1.0;
   int iter = 1;
   bool Lsuccess = false, Usuccess = false;
   
   L = 0.0;
   U = SML::INFTY;
   while (1) {
      //calculate sup_g <g,p> at eta = test_point
      sup_subgrad = 0.0;
      for (vector<vector<Line> >::iterator it = Ls.begin(); it < Ls.end(); it++) {
         double max_slope = 0.0, max_value = -SML::INFTY;
         
         for (vector<Line>::iterator itt = it->begin(); itt < it->end(); itt++) {
            double f = itt->bias + test_point * itt->slope;
            if (f > max_value || (f == max_value && itt->slope > max_slope)) {
               max_slope = itt->slope;
               max_value = f;
               
            }
         }
         sup_subgrad += max_slope;
      }
      sup_subgrad /= numElements;
      sup_subgrad += lambda_wp + test_point * h;
      
      if (sup_subgrad > 0) {
         U = test_point;
         test_point /= 2.0;
         Usuccess = true;
      } else if (sup_subgrad < 0) {
         L = test_point;
         test_point *= 2;
         Lsuccess = true;
      } else {
         L = U = test_point;
         break;
      }
      if (iter > MAX_LU_SEARCH || (Lsuccess && Usuccess)) {
         break;
      }
      iter++;
   }
}

int CMultilabel_Linesearch::FindHinges(double &L, double &U,
                                       bool calculate_value) {
   
   int totalLines = 0;
   for (vector<vector<Line> >::iterator it = Ls.begin(); it < Ls.end(); it++) {
      if (calculate_value) {
         for (vector<Line>::iterator itt = it->begin(); itt < it->end(); itt++) {
            itt->value = itt->bias + L * itt->slope;
         }
      }
      //sort by non-ascending value of each line
      sort(it->begin(), it->end(), descent);
      
      //find max envelope
      deque<lineWithLabel> S;
      lineWithLabel maxLine;
      maxLine.step = L;
      maxLine.slope = (it->front()).slope;
      maxLine.bias = (it->front()).bias;
      //need not to store the dummy label pair (-1, -1)
      if ((it->front()).labels.label != -1 && (it->front()).labels.trueLabel
          != -1) {
         maxLine.labels = vector<LabelPair> (1, (it->front()).labels);
      }
      
      S.push_front(maxLine);
      
      for (vector<Line>::iterator itt = it->begin() + 1; itt < it->end(); itt++) {
         double etaTest = SML::INFTY, eta = 0.0, topSlope = 0.0;
         while (!S.empty()) {
            deque<lineWithLabel>::iterator topLine = S.begin();
            topSlope = topLine->slope;
            eta = topLine->step;
            //handle identical lines
            if (!(abs(itt->slope - topLine->slope) > MYZERO) && 
                !(abs(itt->bias - topLine->bias) > MYZERO)) {
               if (itt->labels.label != -1 && itt->labels.trueLabel != -1) {
                  (topLine->labels).push_back(itt->labels);
               }
               break;
            } else {
               //calculate intersection
               etaTest = (topLine->bias - itt->bias) / (itt->slope - topLine->slope);
               
               if ((etaTest <= eta && etaTest > L) || 
                   (etaTest == L && itt->slope > topLine->slope)) {
                  S.pop_front();
               } else {
                  break;
               }
            }
         }
         if ((etaTest < U && etaTest > L) || (etaTest == L && itt->slope > topSlope)) {
            lineWithLabel maxLine;
            maxLine.step = etaTest;
            maxLine.slope = itt->slope;
            maxLine.bias = itt->bias;
            if (itt->labels.label != -1 && itt->labels.trueLabel != -1) {
               maxLine.labels = vector<LabelPair> (1, itt->labels);
            }
            
            S.push_front(maxLine);
         }
      }
      
      totalLines += S.size();
      if (Ss.size() < uint(numElements)) {
         Ss.push_back(S);
      } else {
         Ss.at(it - Ls.begin()) = S;
      }
      
   }
   return totalLines;
}

int CMultilabel_Linesearch::FindHinges(TheMatrix& f, TheMatrix& df, double L,
                                       double U, bool store_activelabel) {
   
   int totalLines = 0;
   const vector<vector<double> >& Y = _data->labels();
   Line oneLine;
   vector<Line> lines;
   vector<Line>::iterator it_lines;
   //a horizontal line, y = 0. note we set trueLabel=-1
   //and label=-1 since this line is not associated with any label.
   //Line hLine = Line(0.0, LabelPair(-1, -1), 0.0, 0.0);
   bool first = true;
   for (int i = 0; i < numElements; i++) {
      deque<lineWithLabel> S;
      f.CreateMatrixRowView(i, f_i);
      df.CreateMatrixRowView(i, df_i);
      if (!first) {
         //adjust space for the vector lines.
         lines.resize((numClass - Y[i].size()) * Y[i].size() + 1);
         it_lines = lines.begin();
         *(it_lines) = hLine;
         it_lines++;
         
      } else {
         //add a horizontal line for max(0,.) operation of hinge loss
         lines.push_back(hLine);
      }
      //create linear lines, and set 'Line.value' to be the value
      //evaluated at 0
      for (int k = 0; k < numClass; k++) {
         bool isTrueLabel = false;
         //check whether k is in the correct label set
         for (uint j = 0; j < Y[i].size(); j++) {
            if (k + 1 == Y[i][j]) {
               isTrueLabel = true;
               break;
            }
         }
         
         //if k is not in the correct label set
         //then add numCorrectLabel number of linear lines:
         //l := bias + stepsize*slope, where the values of bias and slope
         //depend on the pair of label (k, z), where z
         //is in the correct label set
         if (!isTrueLabel) {
            for (uint j = 0; j < Y[i].size(); j++) {
               double f_k, df_k, f_y, df_y;
               f_i->Get(k, f_k);
               df_i->Get(k, df_k);
               f_i->Get(int(Y[i][j]) - 1, f_y);
               df_i->Get(int(Y[i][j]) - 1, df_y);
               //bias
               oneLine.bias = (1.0 + f_k - f_y);
               //slope
               oneLine.slope = df_k - df_y;
               //value at lower bound
               oneLine.value = oneLine.bias + L * oneLine.slope;
               //associated label pair. note k is shifted by 1
               oneLine.labels = LabelPair(int(Y[i][j]), k + 1);

               //this is done to prevent numerical precision error
               oneLine.value *= (abs(oneLine.value) > MYZERO);
               oneLine.slope *= (abs(oneLine.slope) > MYZERO);
               oneLine.bias *= (abs(oneLine.bias) > MYZERO);
               if (first) {
                  lines.push_back(oneLine);
               } else {
                  *(it_lines) = oneLine;
                  it_lines++;
               }
            }
         }
      }
      
      //sort by non-ascending value of each line
      sort(lines.begin(), lines.end(), descent);
      
      //find max envelope
      lineWithLabel maxLine;
      maxLine.step = L;
      maxLine.slope = (lines.front()).slope;
      maxLine.bias = (lines.front()).bias;
      //need not to store the dummy label pair (-1, -1)
      if (store_activelabel) {
         if ((lines.front()).labels.label != -1
             && (lines.front()).labels.trueLabel != -1) {
            maxLine.labels = vector<LabelPair> (1, (lines.front()).labels);
         }
      }
      S.push_front(maxLine);
      
      for (vector<Line>::iterator it = lines.begin() + 1; it < lines.end(); it++) {
         double etaTest = SML::INFTY, eta = 0.0, topSlope = 0.0;
         while (!S.empty()) {
            deque<lineWithLabel>::iterator topLine = S.begin();
            topSlope = topLine->slope;
            eta = topLine->step;
            //handle identical lines
            if (!(abs(it->slope - topLine->slope) > MYZERO) && !(abs(it->bias - topLine->bias) > MYZERO)) {
               if (store_activelabel) {
                  if (it->labels.label != -1 && it->labels.trueLabel != -1) {
                     (topLine->labels).push_back(it->labels);
                  }
               }
               break;
            } else {
               //calculate intersection
               etaTest = (topLine->bias - it->bias) / (it->slope - topLine->slope);
               if ((etaTest <= eta && etaTest > L) || (etaTest == L && it->slope > topLine->slope)) {
                  S.pop_front();
               } else {
                  break;
               }
            }
            
         }
         if ((etaTest < U && etaTest > L) || (etaTest == L && it->slope > topSlope)) {
            lineWithLabel maxLine;
            maxLine.step = etaTest;
            maxLine.slope = it->slope;
            maxLine.bias = it->bias;
            if (store_activelabel) {
               if (it->labels.label != -1 && it->labels.trueLabel != -1) {
                  maxLine.labels = vector<LabelPair> (1, it->labels);
               }
            }
            S.push_front(maxLine);
         }
      }
      
      totalLines += S.size();
      if (Ss.size() < uint(numElements)) {
         Ss.push_back(S);
      } else {
         Ss.at(i) = S;
      }
      
      if (first) {
         first = false;
      }
      
   } //end of i
   
   return totalLines;
   
}

double CMultilabel_Linesearch::FindOptStep(const double &h,
                                           const double &lambda_wp, int &totalLines) {
   
   return FindOptStep(h, lambda_wp, totalLines, activeLine);
}

double CMultilabel_Linesearch::FindOptStep(const double &h,
                                           const double &lambda_wp, int &totalLines,
                                           vector<lineWithLabel> &activeLine) {
   
   double rho = 0.0, prev_rho = 0.0;
   rho = lambda_wp;
   if (activeLine.size() != uint(numElements)) {
      activeLine.clear();
   }
   for (vector<deque<lineWithLabel> >::iterator it = Ss.begin(); it < Ss.end(); it++) {
      //operate the stack of lines, S, in reverse order
      deque<lineWithLabel>::iterator itt = it->end() - 1;
      //calculating initial slope, i.e., the supremum of subgradients at 0
      rho += itt->slope / numElements;
      //store active lines for later use
      if (activeLine.size() < uint(numElements)) {
         activeLine.push_back(*itt);
      } else {
         activeLine.at(it - Ss.begin()) = *itt;
      }

      //keep track of lines
      if (vec_iterator.size() < uint(numElements)) {
			vec_iterator.push_back(--itt);
      } else {
         vec_iterator.at(it - Ss.begin()) = --itt;
      }
      totalLines--;
   }
   
   CVector_Operations vecOperation(Ss, vec_iterator);
   vector<lineWithLabelAndIndex> nextHinge;
   
   double eta = 0.0;
   // supremum of subgradient at lower bound
   double sup_subgrad = rho;
   int iter = 0;
	while (sup_subgrad < 0) {
       prev_rho = rho;
       if (totalLines == 0) {
          eta = SML::INFTY;
          break;
       }
       
       //find next subdifferentiable point using a priority_queue, cf. basic_utilities.hpp
       nextHinge = vecOperation.Min(Ss, vec_iterator);

       //update supremum of subgradient
       for (vector<lineWithLabelAndIndex>::iterator it = nextHinge.begin(); it
               < nextHinge.end(); it++) {
          //note (it->line) is a iterator of deque<lineWithLabel>
          eta = (it->line)->step;
          //deduct previous slope
          rho -= (it->line + 1)->slope / numElements;
          //add current slope
          rho += (it->line)->slope / numElements;
          totalLines--;
       }
       
		// supremum of subgradients at current eta
       sup_subgrad = rho + eta * h;
       iter++;
       if (sup_subgrad < 0) {
          // update active label information
          for (vector<lineWithLabelAndIndex>::iterator it = nextHinge.begin(); it
                  < nextHinge.end(); it++) {
             activeLine.at(it->index) = *(it->line);
          }
       }
	}
    
	//	printf("ls iter: %d\n", iter);
	double proposedEta = -prev_rho / h;
    
	if (eta > proposedEta) {
       // in this case we are not sitting at any hinge
       eta = proposedEta;
       //printf("exact step: %e\n", eta);
	} else {
       // in this case we are sitting at at least one hinge
       // we save active label of the last active linear line
       // for later use, e.g. calculate the loss
       for (vector<lineWithLabelAndIndex>::iterator it = nextHinge.begin(); it
               < nextHinge.end(); it++) {
          activeLine.at(it->index).slope = (it->line)->slope;
          activeLine.at(it->index).bias = (it->line)->bias;
          
          ((activeLine.at(it->index)).labels).insert(
             (activeLine.at(it->index)).labels.end(),
             (it->line)->labels.begin(), (it->line)->labels.end());
       }
	}
	return eta;
}


double CMultilabel_Linesearch::Linesearch(const TheMatrix &w,
                                          const TheMatrix &p, 
                                          TheMatrix &f, 
                                          TheMatrix &df, 
                                          const double &lambda,
                                          bool calculate_loss, 
                                          bool check_descent) 
{
   double lambda_wp = 0.0;
   double eta = 0.0, sup_subgrad = 0.0;
   p.Dot(w, lambda_wp);
   lambda_wp *= lambda;
   
   if (check_descent) {
      eta = 0.0;
      //first of all, check whether p is a descent direction
      //by looking at the supremum  of the subdifferential
      //ConstructLines(f, df); //saving each line consumes too much memory
      //ComputeLossAndSupGrad(eta, sup_subgrad, error);
      ComputeLossAndSupGrad(eta, sup_subgrad, error, f, df);
      
      sup_subgrad /= numElements;
      sup_subgrad += lambda_wp;
      error /= numElements;
      
      if (sup_subgrad >= 0) {
         return eta;
      }
      
   }
   
   double h = 0.0;
   p.Norm2(h);
   h = h * h * lambda;
   
   int totalLines = FindHinges(f, df);
   
   eta = FindOptStep(h, lambda_wp, totalLines);
   
   // calculate loss if required
   if (calculate_loss) 
   {
      error = 0.0;
      for (vector<lineWithLabel>::iterator it = activeLine.begin(); it
              < activeLine.end(); it++) {
         error += max(0.0, it->bias + eta * it->slope);
      }
      error = error / numElements;
   }
   
   return eta;
}

double CMultilabel_Linesearch::Linesearch(const TheMatrix &p, 
                                          TheMatrix &f,
                                          TheMatrix &df, 
                                          const double &lambda, 
                                          vector<lineWithLabel> &activeLine,
                                          bool calculate_loss) 
{
	int totalLines = FindHinges(f, df);
    
	double h = 0.0, lambda_wp = 0.0;
    
	p.Norm2(h);
	h = h * h * lambda;
	TheMatrix &w = _model->GetW();
	p.Dot(w, lambda_wp);
	lambda_wp *= lambda;
    
	double eta = FindOptStep(h, lambda_wp, totalLines, activeLine);
    
	// calculate loss if required
	if (calculate_loss) 
    {
       error = 0.0;
       for (vector<lineWithLabel>::iterator it = activeLine.begin(); it
               < activeLine.end(); it++) {
          error += max(0.0, it->bias + eta * it->slope);
       }
       error = error / numElements;
	}
    
	return eta;
}

double CMultilabel_Linesearch::LinesearchInTwoSteps(const TheMatrix &p,
                                                    TheMatrix &f, 
                                                    TheMatrix &df, 
                                                    const double &lambda,
                                                    vector<lineWithLabel> &activeLine, 
                                                    bool calculate_loss) 
{
   double L = 0, U = SML::INFTY;
   double h = 0.0, lambda_wp = 0.0;
   
   p.Norm2(h);
   h = h * h * lambda;
   
   TheMatrix &w = _model->GetW();
   p.Dot(w, lambda_wp);
   lambda_wp *= lambda;
   
   ConstructLines(f, df);
   
   FindLAndU(h, lambda_wp, L, U);
   
   int totalLines = FindHinges(L, U);
   
   //	printf("L: %e, U: %e, total: %d \n", L, U, totalLines);
   
   double eta = FindOptStep(h, lambda_wp, totalLines, activeLine);
   
   // calculate loss if required
   if (calculate_loss) 
   {
      error = 0.0;
      for (vector<lineWithLabel>::iterator it = activeLine.begin(); it
              < activeLine.end(); it++) {
         error += max(0.0, it->bias + eta * it->slope);
      }
      error = error / numElements;
   }
   
   return eta;
}

double CMultilabel_Linesearch::InexactLinesearch(const TheMatrix &w,
                                                 const TheMatrix &p, TheMatrix &f, TheMatrix &df, const double &lambda,
                                                 double eta_0, double c, double a) {
   // now check if p is a descent direction by looking at sup_g <g,p>
   double sup_subgrad = 0.0, lambda_wp = 0.0;
   p.Dot(w, lambda_wp);
   lambda_wp *= lambda;
   
   //calculate sup_g <g,p> at eta = 0
   sup_subgrad = 0.0;
   error = 0.0;
   double eta = 0.0;
   ComputeLossAndSupGrad(eta, sup_subgrad, error, f, df);
   
   //	for (vector<deque<lineWithLabel> >::iterator it = Ss.begin(); it < Ss.end(); it++) {
   //		deque<lineWithLabel>::iterator itt = it->end() - 1;
   //		sup_subgrad += itt->slope;
   //		error += max(0.0, itt->bias);
   //	}
   
   sup_subgrad /= numElements;
   sup_subgrad += lambda_wp;
   error /= numElements;
   
   if (sup_subgrad >= 0) {
      //in this case, p is not a descent direction
      //return zero step
      eta = 0.0;
   } else {
      
      //find hinges, and construct sorted stacks
      FindHinges(f, df, 0, SML::INFTY, true);
      // backtracking line search to find a step size that satisfies the sufficient decrease
      // condition. The objective function takes the form of
      // obj = 0.5lambda|w|^2 + lambda*eta*<w,p> + 0.5lambda*eta^2|p|^2 + averaged losses.
      // Note that we can neglect the first term in the objective since it doesn't depend on
      // step size
      double h = 0.0;
      p.Norm2(h);
      h = h * h * lambda / 2.0;
      eta = eta_0;
      int iter = 1;
      bool success = false;
      double loss_0 = error;
      double dummy = 0.0;
      while (1) {
         ComputeLossAndSupGradUsingSs(eta, dummy, error);
         error /= numElements;
         if (error + eta * lambda_wp + eta * eta * h <= loss_0 + c * eta
             * sup_subgrad) {
            success = true;
            break;
         }
         if (iter > MAX_INEXACT_LINESEARCH) {
            break;
         }
         eta *= a;
         iter++;
         
      }
      
      if (!success) {
         eta = 0.0;
         error = loss_0;
      }
   }
   return eta;
}
