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
 * Authors: Qinfeng Shi (qinfeng.shi@anu.edu.au)
 *          Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (14/04/2008) 
 *
 * Last Updated: (28/10/2008))
 */


#ifndef _SMMMULTICLASSLOSS_HPP_
#define _SMMMULTICLASSLOSS_HPP_

#include <string>
#include "sml.hpp"
#include "data.hpp"
#include "seqmulticlassdata.hpp"
#include "loss.hpp"
#include "model.hpp"


/** Class for semi markov model for sequence segmentation and classification
 *  
 *  Reference:
 *  1. Q. Shi , L. Wang, L. Cheng and A.J. Smola, 
 *     "Discriminative Human Action Segmentation and Recognition using Semi-Markov Model", 
 *     In IEEE CVPR 08, To Appear , 2008.
 *  2. Q. Shi, Y. Altun, A.J. Smola, and S.V.N. Vishwanathan,
 *     "Semi-Markov Models for Sequence Segmentation",
 *     in EMNLP, 2007.
 */
class CSMMMulticlassLoss : public CLoss 
{    
   protected:
      /** whether to do crossvalidation
       */
      bool iscrossvalidation;

      TheMatrix *tphi_1;
      TheMatrix *tphi_2;


      /** Pointer to data
       */
      CSeqMulticlassData* _data;
      
      
      /** Number of examples in data
       */
      unsigned int m;    
      
      /** Starting index of phi_1 raw feature, which will be extended to multiclass feature.
       */
      unsigned int startIndexPhi1;
      
      /** Starting index of phi_2 raw feature, which will be extended to multiclass feature.
       */
      unsigned int startIndexPhi2;
      
      /** Starting index of phi_3 raw feature, which will be extended to multiclass feature.
       */
      unsigned int startIndexPhi3;
      
      /** Minimum segment duration
       */
      unsigned int minDuration;
      
      /** Maximum segment duration
       */
      unsigned int maxDuration;
      
      /** Maximum sequence length
       */
      unsigned int maxSequenceLength;
      
      /** Iteration number
       */
      unsigned int iterNum;
      
      
      /** find best label WITH label loss (for training)
       */
      void find_best_label(const std::vector<unsigned int> &y,const std::vector<unsigned int> &ylabel, const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, std::vector<unsigned int> &ybar,std::vector<unsigned int> &ybarlabel, double &marginloss, double &labelloss, unsigned int personid, unsigned int classNum);
      
      /** find best label WITH label loss, WITH a grammer(for training)
       */
      void find_best_label_grammer(const std::vector<unsigned int> &y,const std::vector<unsigned int> &ylabel, const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, std::vector<unsigned int> &ybar,std::vector<unsigned int> &ybarlabel, double &marginloss, double &labelloss, unsigned int personid, unsigned int classNum);
      
      /** find best label WITHOUT label loss (for testing)
       */
      void find_best_label(const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, std::vector<unsigned int> &ybar, std::vector<unsigned int> &ybarlabel, double &marginloss, unsigned int personid, unsigned int classNum);
      
      /** find best label WITHOUT label loss, WITH a grammer(for testing)
       */
      void find_best_label_grammer(const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, std::vector<unsigned int> &ybar, std::vector<unsigned int> &ybarlabel, double &marginloss, unsigned int personid, unsigned int classNum);
      
      /** compute precision and recall
       */
      void PrecRec(const std::vector<unsigned int> &y, const std::vector<unsigned int> &ybar, double &prec, double& rec);
      
      /** "Partial" and decomposed loss proposed in ref[2]
       */
      double PartialDelta(const unsigned int b1, const unsigned int b2, const std::vector<unsigned int> y, const std::vector<unsigned int> &ylabel, unsigned int classid, unsigned int length);
      
      /** "Partial" and decomposed segment-wise loss proposed by ref[2]
       */
      double PartialSegDelta(const unsigned int ybar_i, const std::vector<unsigned int> &y);
      
      /** "Partial" and decomposed segment-wise loss in ref[2]
       */
      double PartialNearestSegDelta(const unsigned int ybar_i, const std::vector<unsigned int> y, unsigned int length);
      
      /** "Partial" and decomposed label point-wise loss proposed by ref[1]
       */
      double PartialPointDelta(const unsigned int b1, const unsigned int b2, std::vector<unsigned int> y, const std::vector<unsigned int> &ylabel, unsigned int classid);
      
      
      
      /** Full label loss proposed by ref[2]
       */
      double Delta(const std::vector<unsigned int> &ybar, const std::vector<unsigned int> &y);
      double NearestSegDelta(const std::vector<unsigned int> &ybar, const std::vector<unsigned int> &y,unsigned int length);
      double AllDelta(const std::vector<unsigned int> &ybar, const std::vector<unsigned int> &y, 
                      const std::vector<unsigned int> &ybarlabel, const std::vector<unsigned int> &ylabel, const int length);
      
      /** Compute decision function values and predict labels
       *
       *  @param model [read] Model object containing a weight vector
       *  @param ybar [write] Predicted labels
       *  @param f_ybar [write] Decision function values
       */
      virtual void DecisionAndPrediction(CModel* model, int* ybar, double* f_ybar);
      
      /**  Covert boundary to sequence 
       */
      std::vector <unsigned int> Boundry2StatSequence(const std::vector<unsigned int> segs,const std::vector<unsigned int> labels,unsigned int length);
      
      /** Get accuracy
       */
      double Accuracy(const std::vector<unsigned int> target,const std::vector<unsigned int> predict);
      double Labelloss(const std::vector<unsigned int> target,const std::vector<unsigned int> predict);
      
      
      
   public:
      CSMMMulticlassLoss() {}
      CSMMMulticlassLoss(CModel* &model, CSeqMulticlassData* &data);
      virtual ~CSMMMulticlassLoss();
      
      // Methods
      virtual void ComputeLoss(double& loss);
      virtual void ComputeLoss(std::vector<unsigned int> y, std::vector<unsigned int> ylabel, std::vector<unsigned int> ybar, std::vector<unsigned int> ybarlabel, const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, double & marginloss, double & labelloss, int flag);
      virtual void ComputeLossAndGradient(double& loss, TheMatrix& grad);
      virtual void Evaluate(CModel* model);
      virtual void EvaluateAll(CModel* model,double & weighted_prec,double & weighted_rec, double& weighted_avgf1,
                               std::vector <unsigned int> &sseqStdAll,std::vector <unsigned int> &sseqBarAll,std::vector<unsigned int> &pid,
                               std::vector<unsigned int> & pid_bar,double & avgpid, double & voted_acc, double & weighted_acc);
      virtual void EvaluateMark(CModel* model, int mark, double & weighted_prec,double & weighted_rec, double& weighted_avgf1,
                                std::vector <unsigned int> &sseqStdAll,std::vector <unsigned int> &sseqBarAll,std::vector<unsigned int> &pid,
                                std::vector<unsigned int> & pid_bar,double & avgpid, double & voted_acc, double & weighted_acc);
      
};

#endif
