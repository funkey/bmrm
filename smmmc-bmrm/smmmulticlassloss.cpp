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
 * Last Updated: (28/10/2008)
 */


#ifndef _SMMMULTICLASSLOSS_CPP_
#define _SMMMULTICLASSLOSS_CPP_

#define MAX(x,y)      ((x) < (y) ? (y) : (x))
#define MIN(x,y)      ((x) > (y) ? (y) : (x))

//const double yeta = 0.01;
//const double yeta = 0.1;
const double yeta = 1;
//int lastDuration = 30; % for synthetic set
int lastDuration = -1;
bool is_first_phi1_used = true;
bool is_single_action_persequence = false;// use this as a grammer
bool is_pseudo_sgd = false; //In every iteration,the gradient and loss are only caculated on a random subset of the training examples
unsigned int exam_per_iter = 1;
//seg, label, seg-nearest
double lossw[3] = {0.5,1,0.5} ;


#include "common.hpp"  // WriteFile()
#include "smmmulticlassloss.hpp"
#include "configuration.hpp"
#include <sstream>

using namespace std;


double acc_all,acc_train,acc_test,acc_validation;
double recall_all,precision_all,f1score_all;
double recall_train,precision_train,f1score_train;
double recall_test,precision_test,f1score_test;
double recall_valid,precision_valid,f1score_valid;
double acc_all_vote,acc_train_vote,acc_test_vote,acc_validation_vote;
vector <unsigned int> allStd_test,allBar_test,allStd_valid,allBar_valid,allStd_train,allBar_train,allStd,allBar;
vector <unsigned int> allStd_pid,allPred_pid,trainStd_pid,trainPred_pid,testStd_pid,testPred_pid,validStd_pid,validPred_pid;
double acc_all_pid,acc_train_pid,acc_test_pid,acc_validation_pid;


/* find a certain value in a vector */
int find(unsigned int val,vector <unsigned int> &a)
{
    unsigned int num = (unsigned int)a.size();
    for(unsigned int i=0;i<num;i++)
    {
	if(val==a[i])
	    return i;
    }
    return -1;

}
/* voting */
unsigned int vote(vector <unsigned int> &a)
{
    unsigned int num = (unsigned int)a.size();
    vector <unsigned int> b;
    vector <unsigned int> freq;
    
    for(unsigned int i=0;i<num;i++)
    {
	int pos = find(a[i],b);
	if (pos == -1)
	{
	    b.push_back(a[i]);
	    freq.push_back(1);
	}
	else
	{
	    freq[pos] += 1;
	}
    }
    unsigned int size_freq = freq.size();
    unsigned int majorPos = 0;
    for(unsigned int i=0;i<size_freq;i++)
    {		
	if(freq[i]>freq[majorPos])
	{
	    majorPos = i;
	}
    }
    
    return b[majorPos];
    
}



/** Constructor
 *
 *  @param w    [r/w] external pointer to weight vector
 *  @param data [r/w] pointer to data container
 *  @param additiveLabelLoss [read] flag to determine whether to "add" label loss to loss function or not (default: true)
 */
CSMMMulticlassLoss::CSMMMulticlassLoss(CModel* &model, CSeqMulticlassData* &data)
   : iscrossvalidation(false),
     tphi_1(0),
     tphi_2(0),
     _data(data),
     m(0)          
{
    iterNum = 0;
    maxSequenceLength = 0;
    // get/set problem parameters
    Configuration &config = Configuration::GetInstance();   
    
    if(config.IsSet("SMM_MULTICLASS.minDuration"))
       minDuration = config.GetInt("SMM_MULTICLASS.minDuration");
    else
       minDuration = _data->MinDuration();
    
    if(config.IsSet("SMM_MULTICLASS.maxDuration"))
       maxDuration = config.GetInt("SMM_MULTICLASS.maxDuration");
    else
       maxDuration = _data->MaxDuration();  
    assert(minDuration <= maxDuration);
    
    if(config.IsSet("SMM_MULTICLASS.doCrossValidation"))
       iscrossvalidation = config.GetBool("SMM_MULTICLASS.doCrossValidation");

    m = _data->slice_size();        
    assert(m > 0);
    scalingFactor = 1.0/m;
    maxSequenceLength = _data->MaxSequenceLength(); 
    
    startIndexPhi1 = _data->StartIndexPhi1();
    startIndexPhi2 = _data->StartIndexPhi2();
    startIndexPhi3 = _data->StartIndexPhi3();
    
    // extend raw features to tensor features
    _data->ExtendFeatures();
    
    // create temporary tensor features
    if(_data->tdim() <= 0)
       throw CBMRMException("_data->tdim() must be > 0", "CSMMMulticlassLoss::CSMMMulticlassLoss()");
    tphi_1 = new TheMatrix(1,_data->tdim(), SML::DENSE);
    tphi_2 = new TheMatrix(1,_data->tdim(), SML::DENSE);
    
    // initialize model
    if(not model->IsInitialized())
    {
       model->Initialize(1, _data->tdim(), _data->bias());
       printf("rawfeatureDim:%d\n tensorfeatureDim:%d\n",_data->dim(),_data->tdim());
    }
    
    // keep a pointer to model
    _model = model;
}


CSMMMulticlassLoss::~CSMMMulticlassLoss()
{
    if(tphi_1) delete tphi_1;
    if(tphi_2) delete tphi_2;
}



/**  Compute loss
 */
void CSMMMulticlassLoss::ComputeLoss(double& loss)
{
    TheMatrix dummygrad(_model->GetNumOfW(), _model->GetDimOfW(), SML::DENSE);
    ComputeLossAndGradient(loss,dummygrad); 
}

/** Flag = 0: marginloss, no label loss. The label loss will always be zero
           1: marginloss, and label loss.
*/
void CSMMMulticlassLoss::ComputeLoss(vector<unsigned int> y, vector<unsigned int> ylabel, vector<unsigned int> ybar, vector<unsigned int> ybarlabel, const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, double & marginloss, double & labelloss, int flag)
{
    unsigned int i;
    double w_dot_phi1 = 0;
    double w_dot_phi2 = 0;
    marginloss = 0;

    unsigned int start;
    if(is_first_phi1_used)
	start = 0;
    else
	start = 1;
    for(i=start; i < ybar.size(); i++)
    {
       _data->TensorPhi1(x.phi_1[ybar[i]],ybarlabel[i],0,tphi_1);
       //tphi_1->Print();
       w.Dot(*(tphi_1), w_dot_phi1);
       marginloss += w_dot_phi1;
       //printf("%d(%d):%2.4f\t",ybar[i],ybarlabel[i],marginloss);
    }	
    for(i=1;i<ybar.size();i++)
    {
       int vb = 0;
       _data->TensorPhi2(x.phi_2[ybar[i-1]][ybar[i]-ybar[i-1]-1], ybarlabel[i-1], ybarlabel[i], 0,vb,tphi_2);
       w.Dot(*(tphi_2), w_dot_phi2);
       marginloss += w_dot_phi2;
    }
    
    if(ybar.size() > 0)
    {
       
       //grad.Add(*(X[i].phi_2[ybar[ybar.size()-1]][X[i].len-1 - ybar[ybar.size()-1]-1]));////       
       _data->TensorPhi2(x.phi_2[ybar[ybar.size()-1]][x.len - ybar[ybar.size()-1]-1 ], ybarlabel[ybar.size()-1], 0, 0,0,tphi_2);
       w.Dot(*(tphi_2), w_dot_phi2);
       marginloss += w_dot_phi2;
    }
    
    //vector <unsigned int> yss = Boundry2StatSequence(y,ylabel,x.len);
    //vector <unsigned int> ybarss = Boundry2StatSequence(ybar,ybarlabel,x.len);
    //labelloss = Labelloss(yss,ybarss);
    labelloss = AllDelta(ybar,y,ybarlabel,ylabel,x.len);
}



/**   Compute loss and gradient
 */
void CSMMMulticlassLoss::ComputeLossAndGradient(double& loss, TheMatrix& grad)
{
   iterNum ++;
   TheMatrix &w = _model->GetW();
   loss = 0;
   grad.Zero();
   TheMatrix g(grad, SML::DENSE);
   
   const vector<CSeqMulticlassLabel::seqlabel_struct> &Y = _data->labels();
   const vector<CSeqMulticlassFeature::seqfeature_struct> &X = _data->features();
   
   unsigned int trainExNum = 0;
   vector <int > cvmark = _data->Getcvmark();	
   for(unsigned int i=0; i < m; i++)
   {
      if(cvmark.size()!=0)			
      {
         if(cvmark[i]!=SMM::TRAIN_DATA)
            continue;
      }
      trainExNum ++;
      
      //if(cvmark)
      vector<unsigned int> ybar(X[i].len,0);
      vector<unsigned int> ybarlabel(X[i].len,0);
      double labelloss = 0;
      double marginloss = 0;
      double w_dot_g = 0.0;;
      
      // find best label y' and return the score wrt to y'
      if(verbosity>=2)
      {
         cout <<"ex:"<< i<< endl;fflush(stdout);
      }
      
      if(is_single_action_persequence)
         find_best_label_grammer(Y[i].pos,Y[i].type, X[i], w, ybar, ybarlabel, marginloss, labelloss, 0, _data->getNumOfClass());
      else
         find_best_label(Y[i].pos,Y[i].type, X[i], w, ybar, ybarlabel, marginloss, labelloss, 0, _data->getNumOfClass());
      
      double labelloss_y = 0;
      double marginloss_y = 0;
      double labelloss_ybar = 0;
      double marginloss_ybar = 0;
      
      
      ComputeLoss(Y[i].pos,Y[i].type,ybar,ybarlabel,X[i],w,marginloss_ybar,labelloss_ybar,1);
      if(lossw[0]!=0)
         labelloss+=lossw[0];
      
      if(lastDuration>0)
      {
         marginloss = marginloss_ybar;
         labelloss = labelloss_ybar;
      }
      if(verbosity>=3)
      {					
         ComputeLoss(Y[i].pos,Y[i].type,Y[i].pos,Y[i].type,X[i],w,marginloss_y,labelloss_y,1);
         printf("dp------marginloss:%2.4f---labelloss:%2.4f------\n",marginloss,labelloss);	
         printf("ybar----marginloss:%2.4f---labelloss:%2.4f------\n",marginloss_ybar,labelloss_ybar);
         printf("y-------marginloss:%2.4f---labelloss:%2.4f------\n",marginloss_y,labelloss_y);			
         if(abs(labelloss_ybar-labelloss)>1e-5)
         {
            printf("labelloss doesn't match!\n");
            //exit(0);
         }
         if(abs(marginloss_ybar-marginloss)>1e-5)
         {
            printf("marginloss_ybar_dp:%2.4f != marginloss_ybar_computeLoss:%2.4f\n",marginloss,marginloss_ybar);
            printf("marginloss doesn't match!\n");
         }
      }
      
      // construct the gradient vector for the part of true y
      const vector<unsigned int> &y = Y[i].pos;
      const vector<unsigned int> &ylabel = Y[i].type;
      g.Zero();
      
      for(unsigned int j=0; j < y.size(); j++)
      {
         //g.Add(*(X[i].phi_1[y[j]]));
         //g.Add(*(X[i].phi_2[y[j-1]][y[j]-y[j-1]-1]));
         _data->TensorPhi1(X[i].phi_1[y[j]],ylabel[j],0,tphi_1);
         g.Add(*tphi_1);
         if(j > 0)
         {
            _data->TensorPhi2(X[i].phi_2[y[j-1]][y[j]-y[j-1]-1], ylabel[j-1], ylabel[j], 0,0,tphi_2);
            g.Add(*tphi_2);			
         }
      }
      if(y.size() > 0)
      {
         //g.Add(*(X[i].phi_2[y[y.size()-1]][X[i].len-1 - y[y.size()-1]-1]));////
         _data->TensorPhi2(X[i].phi_2[y[y.size()-1]][X[i].len - y[y.size()-1]-1 ], ylabel[y.size()-1], 0,0,0,tphi_2);
         g.Add(*tphi_2);
      }
      
      // for predicted y'
      for(unsigned int j=0; j < ybar.size(); j++)
      {  
         //grad.Add(*(X[i].phi_1[ybar[j]]));                         
         //grad.Add(*(X[i].phi_2[ybar[j-1]][ybar[j]-ybar[j-1]-1]));
         _data->TensorPhi1(X[i].phi_1[ybar[j]],ybarlabel[j],0,tphi_1);
         grad.Add(*tphi_1);
         if(j>0)			
         {
            _data->TensorPhi2(X[i].phi_2[ybar[j-1]][ybar[j]-ybar[j-1]-1], ybarlabel[j-1], ybarlabel[j], 0,0,tphi_2);
            grad.Add(*tphi_2); ////			
         }
      }
      if(ybar.size() > 0)
      {
         //grad.Add(*(X[i].phi_2[ybar[ybar.size()-1]][X[i].len-1 - ybar[ybar.size()-1]-1]));
         _data->TensorPhi2(X[i].phi_2[ybar[ybar.size()-1]][X[i].len - ybar[ybar.size()-1]-1 ], ybarlabel[ybar.size()-1], 0, 0,0,tphi_2);
         grad.Add(*tphi_2);
      }
      grad.Minus(g);
      
      
      // accumulate the loss
      w.Dot(g, w_dot_g);	
      loss = loss - w_dot_g + marginloss + labelloss;    
      
   }
   scalingFactor = 1.0/trainExNum;
   grad.Scale(scalingFactor);	
   loss *= scalingFactor;        
   
   if(verbosity)
   {
      double gnorm = 0.0;
      grad.Norm2(gnorm);
      cout << "gradient norm=" << gnorm << endl;
   }
   //Evaluate(_model);
}


/** find best label (with label loss): g(w) := max_y' <w,\phi(x,y')> + Delta(y', y)
 *
 *  @param x [read] sequence
 *  @param y [read] actual label for x
 *  @param w [read] weight vector
 *  @param ybar [write] found best label
 *  @param marginloss [write] margin loss <w,\Phi(x,y')> w.r.t to best y'
 *  @param labelloss [write] label loss \Delta(y',y) w.r.t. to best y'
 *
 */
void CSMMMulticlassLoss::find_best_label(const vector<unsigned int> &y,const vector<unsigned int> &ylabel, const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar,vector<unsigned int> &ybarlabel, double &marginloss, double &labelloss, unsigned int personid, unsigned int classNum)
{
    // reset return values
    marginloss = 0;
    labelloss = 0;
    ybar.clear();
    ybarlabel.clear();
    
    /** The margin value vector used in dynamic programming
     */
    vector< vector<double> > M (x.len+1,vector<double> (classNum,0));
    
    /** The label loss value vector used in dynamic programming
     */
    vector< vector<double> > L (x.len+1,vector<double> (classNum,0));
    
    /** The back pointers vector used in dynamic programming to retrieve the optimal path
     */
    // The positions
    vector< vector<int> > A (x.len+1,vector<int> (classNum,-1));
    // The class labels
    vector< vector<int> > C (x.len+1,vector<int> (classNum,0));
    
    
    double maxval = -SML::INFTY;
    double w_dot_phi1 = 0;
    double w_dot_phi2 = 0;
    double marginval = 0;
    double labelval = 0;
    unsigned int right = 0;
    unsigned int left = 0;
    unsigned int start = 0;
    unsigned int end = 0;
    unsigned int classID = 0;
    unsigned int classIDPrev = 0;
    
    double sum = 0;
    
    // compute DP statistics for positions 1 to len-1
    //L[0] += y.size()-2;
    //A[1] = 0;
    for(classID=0;classID<classNum;classID++)
    {
	A[1][classID] = 0;
	//C[1][classID] = 0;
    }
    
    //debug
    
    
    //printf("x.len:%d",x.len);
    if(is_first_phi1_used)
    {
	right =0;
	for(classID=0; classID < classNum; classID++)
	{
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    _data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
	    //tphi_1->Print();
	    w.Dot(*(tphi_1), w_dot_phi1);
	    marginval = w_dot_phi1;   					
	    sum = marginval;
	    if(sum > maxval)
	    {
		M[right][classID] = marginval;			
		maxval = sum;
	    }
	}
    }
    
    for(right=1; right < x.len+1; right++)        
    {
	
	for(classID=0; classID < classNum; classID++)
	{		
	    // \Phi = (phi1, phi2[left,right])
	    // <w, \Phi> = <w,phi1> + <w,phi[left,right]>                
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    
	    //w.Dot(*(x.phi_1[right]), w_dot_phi1);
	    //printf("pos:%d,classid:%d ",right,classID);fflush(stdout);
	    //x.phi_1[right]->Print();
	    if(right < x.len)
	    {
		_data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
		//tphi_1->Print();
		w.Dot(*(tphi_1), w_dot_phi1);
	    }		
	    
	    start = max(0,int(right-maxDuration));
	    //end = right;//-minDuration+1;
	    
	    if(lastDuration>0)
	    {			
		unsigned int lastpos = x.len-lastDuration+1 ;
		end = MIN(right,lastpos);
	    }
	    else
		end = right;

	    for(left=start; left < end; left++)
	    {
		//labelval = PartialDelta(right,y);   
		
		//printf(" b1:%d,b2:%d (%d) loss:%f ",left,right,classID,labelval);
		//if(right>=61)
		//{
		//    getchar();					
		//}
		for(classIDPrev=0;classIDPrev<classNum;classIDPrev++)
		{
		    labelval = PartialDelta(left,right,y,ylabel,classIDPrev,x.len);
		    assert( (labelval<=x.len) && (labelval>=0) );
		    int vb = 0;
		    _data->TensorPhi2(x.phi_2[left][right-left-1], classIDPrev, classID, 0,vb,tphi_2);
		    w.Dot(*(tphi_2), w_dot_phi2); 
		    marginval = w_dot_phi1 + w_dot_phi2;   
		    sum = M[left][classIDPrev]+marginval + L[left][classIDPrev]+labelval;
		    if(sum > maxval)
		    {
			A[right][classID] = left;
			C[right][classID] = classIDPrev;
			M[right][classID] = M[left][classIDPrev] + marginval;
			L[right][classID] = L[left][classIDPrev] + labelval;
			maxval = sum;
		    }
		}
	    }
	}
    }

        
    // get optimal path (i.e. segmentation)        
    unsigned int pos,prepos,classid,preclassid;
    pos = A[x.len][0];
    classid = C[x.len][0];
    
    if(lastDuration>0)
    {	
	pos = x.len-lastDuration;
	classid = 0;
    }
    ybar.push_back(pos);
    ybarlabel.push_back(classid);
    
    prepos = pos;
    preclassid = classid;
	
    while(A[pos][classid] >= 0)
    {                
	pos = A[prepos][preclassid];
	classid = C[prepos][preclassid];
	ybar.push_back(pos);//positions
	ybarlabel.push_back(classid);//class labels		
	
	//printf("%d(%d):%2.4f ",pos,classid,L[pos][classid]);fflush(stdout);
	prepos = pos;
	preclassid = classid;
    }
    
    marginloss = M[x.len][0];
    labelloss = L[x.len][0];
    //printf("finished back track\n labelloss:%3.4f,marginloss:%3.4f\n",labelloss,marginloss);fflush(stdout);
    reverse(ybar.begin(), ybar.end());
    reverse(ybarlabel.begin(), ybarlabel.end());
    
    //printf("reversed\n");fflush(stdout);
    unsigned int i;
    if(verbosity>=2)
    {
	printf("y:   ");
	for(i=0;i<y.size();i++)
	{
	    printf("%d(%d) ",y[i],ylabel[i]);
	}
	fflush(stdout);
	printf("\nybar:");
	for(i=0;i<ybar.size();i++)
	{
	    printf("%d(%d) ",ybar[i],ybarlabel[i]);
	}
	fflush(stdout);
	printf("\nmargin:%f, loss:%f, totalloss:%f\n",marginloss,labelloss,marginloss+labelloss);
    }
}

/** find best label with a grammer(with label loss): g(w) := max_y' <w,\phi(x,y')> + Delta(y', y)
 *
 *  @param x [read] sequence
 *  @param y [read] actual label for x
 *  @param w [read] weight vector
 *  @param ybar [write] found best label
 *  @param marginloss [write] margin loss <w,\Phi(x,y')> w.r.t to best y'
 *  @param labelloss [write] label loss \Delta(y',y) w.r.t. to best y'
 *
 */
void CSMMMulticlassLoss::find_best_label_grammer(const vector<unsigned int> &y,const vector<unsigned int> &ylabel, const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar,vector<unsigned int> &ybarlabel, double &marginloss, double &labelloss, unsigned int personid, unsigned int classNum)
{
    // reset return values
    marginloss = 0;
    labelloss = 0;
    ybar.clear();
    ybarlabel.clear();
    
    /** The margin value vector used in dynamic programming
     */
    
    vector< vector<double> > M (x.len+1,vector<double> (classNum,0));
    
    /** The label loss value vector used in dynamic programming
     */
    vector< vector<double> > L (x.len+1,vector<double> (classNum,0));
    
    /** The back pointers vector used in dynamic programming to retrieve the optimal path
     */
    // The positions
    vector< vector<int> > A (x.len+1,vector<int> (classNum,-1));
    // The class labels
    vector< vector<int> > C (x.len+1,vector<int> (classNum,0));
    
    
    double maxval = -SML::INFTY;
    double w_dot_phi1 = 0;
    double w_dot_phi2 = 0;
    double marginval = 0;
    double labelval = 0;
    unsigned int right = 0;
    unsigned int left = 0;
    unsigned int start = 0;
    unsigned int end = 0;
    unsigned int classID = 0;
    unsigned int classIDPrev = 0;
    
    double sum = 0;
    
    // compute DP statistics for positions 1 to len-1
//         L[0] += y.size()-2;
//         A[1] = 0;
    for(classID=0;classID<classNum;classID++)
    {
	A[1][classID] = 0;
	//C[1][classID] = 0;
    }
    
    //debug
    
    
    //printf("x.len:%d",x.len);
    if(is_first_phi1_used)
    {
	right =0;
	for(classID=0;classID<classNum;classID++)
	{
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    _data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
	    //tphi_1->Print();
	    w.Dot(*(tphi_1), w_dot_phi1);
	    marginval = w_dot_phi1;   					
	    sum = marginval;
	    if(sum > maxval)
	    {
		M[right][classID] = marginval;			
		maxval = sum;
	    }
	}
    }
    
    for(right=1; right < x.len+1; right++)        
    {
	
	for(classID=0;classID<classNum;classID++)
	{		
	    // \Phi = (phi1, phi2[left,right])
	    // <w, \Phi> = <w,phi1> + <w,phi[left,right]>                
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    
	    //w.Dot(*(x.phi_1[right]), w_dot_phi1);
	    //printf("pos:%d,classid:%d ",right,classID);fflush(stdout);
	    //x.phi_1[right]->Print();
	    if(right<x.len)
	    {
		_data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
		//tphi_1->Print();
		w.Dot(*(tphi_1), w_dot_phi1);
	    }		
	    
	    start = max(0,int(right-maxDuration));
	    //end = right;//-minDuration+1;
	    
	    if(lastDuration>0)
	    {			
		unsigned int lastpos = x.len-lastDuration+1 ;
		end = MIN(right,lastpos);
	    }
	    else
		end = right;
	    for(left=start; left < end; left++)
	    {
		classIDPrev = classID;
		labelval = PartialDelta(left,right,y,ylabel,classIDPrev,x.len);
		assert( (labelval<=x.len) && (labelval>=0) );
		int vb = 0;
		_data->TensorPhi2(x.phi_2[left][right-left-1], classIDPrev, classID, 0,vb,tphi_2);
		w.Dot(*(tphi_2), w_dot_phi2); 
		marginval = w_dot_phi1 + w_dot_phi2;   
		sum = M[left][classIDPrev]+marginval + L[left][classIDPrev]+labelval;
		if(sum > maxval)
		{
		    A[right][classID] = left;
		    C[right][classID] = classIDPrev;
		    M[right][classID] = M[left][classIDPrev] + marginval;
		    L[right][classID] = L[left][classIDPrev] + labelval;
		    maxval = sum;
		}
		
		
	    }
	    
	}
    }
    
    // get optimal path (i.e. segmentation)        
    unsigned int pos,prepos,classid,preclassid;
    
    int maxclassid = 0;
    maxval = -SML::INFTY;
    for(unsigned int i=0;i<classNum;i++)
    {
	sum = M[x.len][i] + L[x.len][i];
	if(sum>maxval)
	{
	    maxval = sum;
	    maxclassid = i;
	}
    }
    
    pos = A[x.len][maxclassid];
    classid = C[x.len][maxclassid];
    
    if(lastDuration>0)
    {	
	pos = x.len-lastDuration;
	classid = 0;
    }
    ybar.push_back(pos);
    ybarlabel.push_back(classid);

    prepos = pos;
    preclassid = classid;
    
    
    while(A[pos][classid] >= 0)
    {                
	pos = A[prepos][preclassid];
	classid = C[prepos][preclassid];
	ybar.push_back(pos);//positions
	ybarlabel.push_back(classid);//class labels		
	
	//printf("%d(%d):%2.4f ",pos,classid,L[pos][classid]);fflush(stdout);
	prepos = pos;
	preclassid = classid;
    }
    
    
    marginloss = M[x.len][maxclassid];
    labelloss = L[x.len][maxclassid];
    //printf("finished back track\n labelloss:%3.4f,marginloss:%3.4f\n",labelloss,marginloss);fflush(stdout);
    reverse(ybar.begin(), ybar.end());
    reverse(ybarlabel.begin(), ybarlabel.end());
    
    //printf("reversed\n");fflush(stdout);
    unsigned int i;
    if(verbosity>=2)
    {
	printf("y:   ");
	for(i=0;i<y.size();i++)
	{
	    printf("%d(%d) ",y[i],ylabel[i]);
	}
	fflush(stdout);
	printf("\nybar:");
	for(i=0;i<ybar.size();i++)
	{
	    printf("%d(%d) ",ybar[i],ybarlabel[i]);
	}
	fflush(stdout);
	printf("\nmargin:%f, loss:%f, totalloss:%f\n",marginloss,labelloss,marginloss+labelloss);
    }
}


/** find best label (without label loss): g(w) := max_y' <w,\phi(x,y')>
 *
 *  @param x [read] sequence
 *  @param w [read] weight vector
 *  @param ybar [write] found best label
 *  @param marginloss [write] margin loss <w,\Phi(x,y')> w.r.t to best y'
 */
void CSMMMulticlassLoss::find_best_label(const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar, vector<unsigned int> &ybarlabel, double &marginloss, unsigned int personid, unsigned int classNum)
{
    using namespace std;
    
    // reset return values
    marginloss = 0;        
    ybar.clear();
    ybarlabel.clear();
    
    /** The margin value vector used in dynamic programming
     */	
    vector< vector<double> > M (x.len+1,vector<double> (classNum,0));
    
    /** The back pointers vector used in dynamic programming to retrieve the optimal path
     */
    // The positions
    vector< vector<int> > A (x.len+1,vector<int> (classNum,-1));
    // The class labels
    vector< vector<int> > C (x.len+1,vector<int> (classNum,0));
    
    
    double maxval = -SML::INFTY;
    double w_dot_phi1 = 0;
    double w_dot_phi2 = 0;
    double marginval = 0;
    unsigned int right = 0;
    unsigned int left = 0;
    unsigned int start = 0;
    unsigned int end = 0;
    unsigned int classID = 0;
    unsigned int classIDPrev = 0;
    
    double sum = 0;
    
    // compute DP statistics for positions 1 to len-1
    for(classID=0;classID<classNum;classID++)
    {
	A[1][classID] = 0;
	//C[1][classID] = 0;
    }
    
    if(is_first_phi1_used)
    {
	right =0;
	for(classID=0;classID<classNum;classID++)
	{
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    _data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
	    //tphi_1->Print();
	    w.Dot(*(tphi_1), w_dot_phi1);
	    marginval = w_dot_phi1;   					
	    sum = marginval;
	    if(sum > maxval)
	    {
		M[right][classID] = marginval;			
		maxval = sum;
	    }
	}
    }
    for(right=1; right < x.len+1; right++)
    {
	for(classID=0;classID<classNum;classID++)
	{		
	    // \Phi = (phi1, phi2[left,right])
	    // <w, \Phi> = <w,phi1> + <w,phi[left,right]>                
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    
	    //w.Dot(*(x.phi_1[right]), w_dot_phi1);
	    //printf("pos:%d ",right);
	    //x.phi_1[right]->Print();
	    if(right<x.len)
	    {
		_data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
		w.Dot(*(tphi_1), w_dot_phi1);
	    }		
	    
	    start = max(0,int(right-maxDuration));
	    //end = right;//-minDuration+1;
	    if(lastDuration>0)
	    {			
		unsigned int lastpos = x.len-lastDuration+1 ;
		end = MIN(right,lastpos);
	    }
	    else
		end = right;
	    for(left=start; left < end; left++)
	    {
		for(classIDPrev=0;classIDPrev<classNum;classIDPrev++)
		{					
		    int vb = 0;
		    
		    _data->TensorPhi2(x.phi_2[left][right-left-1], classIDPrev, classID, 0,vb,tphi_2);
		    w.Dot(*(tphi_2), w_dot_phi2); 
		    marginval = w_dot_phi1 + w_dot_phi2;   
		    sum = M[left][classIDPrev]+marginval;
		    if(sum > maxval)
		    {
			A[right][classID] = left;
			C[right][classID] = classIDPrev;
			M[right][classID] = M[left][classIDPrev] + marginval;						
			maxval = sum;
		    }
		}
	    }	    
	}
    }
	
    // get optimal path (i.e. segmentation)        
    unsigned int pos,prepos,classid,preclassid;
    pos = A[x.len][0];
    classid = C[x.len][0];
    if(lastDuration>0)
    {	
	pos = x.len-lastDuration;
	classid = 0;
    }
    ybar.push_back(pos);
    ybarlabel.push_back(classid);
    
    prepos = pos;
    preclassid = classid;
    while(A[pos][classid] >= 0)
    {                
	pos = A[prepos][preclassid];
	classid = C[prepos][preclassid];
	ybar.push_back(pos);//positions
	ybarlabel.push_back(classid);//class labels		

	prepos = pos;
	preclassid = classid;
    }
    marginloss = M[x.len][0];
    reverse(ybar.begin(), ybar.end());
    reverse(ybarlabel.begin(), ybarlabel.end());
}


/** find best label (without label loss): g(w) := max_y' <w,\phi(x,y')>
 *
 *  @param x [read] sequence
 *  @param w [read] weight vector
 *  @param ybar [write] found best label
 *  @param marginloss [write] margin loss <w,\Phi(x,y')> w.r.t to best y'
 */
void CSMMMulticlassLoss::find_best_label_grammer(const CSeqMulticlassFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar, vector<unsigned int> &ybarlabel, double &marginloss, unsigned int personid, unsigned int classNum)
{
    using namespace std;
    
    // reset return values
    marginloss = 0;        
    ybar.clear();
    ybarlabel.clear();
    
    /** The margin value vector used in dynamic programming
     */	
    vector< vector<double> > M (x.len+1,vector<double> (classNum,0));
    
    /** The back pointers vector used in dynamic programming to retrieve the optimal path
     */
    // The positions
    vector< vector<int> > A (x.len+1,vector<int> (classNum,-1));
    // The class labels
    vector< vector<int> > C (x.len+1,vector<int> (classNum,0));
    
    
    double maxval = -SML::INFTY;
    double w_dot_phi1 = 0;
    double w_dot_phi2 = 0;
    double marginval = 0;
    unsigned int right = 0;
    unsigned int left = 0;
    unsigned int start = 0;
    unsigned int end = 0;
    unsigned int classID = 0;
    unsigned int classIDPrev = 0;
    
    double sum = 0;
    
    // compute DP statistics for positions 1 to len-1
    for(classID=0;classID<classNum;classID++)
    {
	A[1][classID] = 0;
	//C[1][classID] = 0;
    }
    
    if(is_first_phi1_used)
    {
	right =0;
	for(classID=0;classID<classNum;classID++)
	{
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    _data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
	    //tphi_1->Print();
	    w.Dot(*(tphi_1), w_dot_phi1);
	    marginval = w_dot_phi1;   					
	    sum = marginval;
	    if(sum > maxval)
	    {
		M[right][classID] = marginval;			
		maxval = sum;
	    }
	}
    }
    for(right=1; right < x.len+1; right++)
    {
	for(classID=0;classID<classNum;classID++)
	{		
	    // \Phi = (phi1, phi2[left,right])
	    // <w, \Phi> = <w,phi1> + <w,phi[left,right]>                
	    maxval = -SML::INFTY;
	    w_dot_phi1 = 0.0;
	    
	    if(right<x.len)
	    {
		_data->TensorPhi1(x.phi_1[right],classID,0,tphi_1);
		w.Dot(*(tphi_1), w_dot_phi1);
	    }		

	    start = max(0,int(right-maxDuration));
	    //end = right;//-minDuration+1;
	    if(lastDuration>0)
	    {			
		unsigned int lastpos = x.len-lastDuration+1 ;
		end = MIN(right,lastpos);
	    }
	    else
		end = right;
	    for(left=start; left < end; left++)
	    {
		classIDPrev = classID;
		int vb = 0;
		_data->TensorPhi2(x.phi_2[left][right-left-1], classIDPrev, classID, 0,vb,tphi_2);
		w.Dot(*(tphi_2), w_dot_phi2); 
		marginval = w_dot_phi1 + w_dot_phi2;   
		sum = M[left][classIDPrev]+marginval;
		if(sum > maxval)
		{
		    A[right][classID] = left;
		    C[right][classID] = classIDPrev;
		    M[right][classID] = M[left][classIDPrev] + marginval;						
		    maxval = sum;
		}
	    }
	    
	}
    }
        
    // get optimal path (i.e. segmentation)        
    unsigned int pos,prepos,classid,preclassid;
    int maxclassid = 0;
    maxval = -SML::INFTY;
    for(unsigned int i=0;i<classNum;i++)
    {
	sum = M[x.len][i];
	if(sum>maxval)
	{
	    maxval = sum;
	    maxclassid = i;
	}
    }
    
    pos = A[x.len][maxclassid];
    classid = C[x.len][maxclassid];
    
    if(lastDuration>0)
    {	
	pos = x.len-lastDuration;
	classid = 0;
    }
    ybar.push_back(pos);
    ybarlabel.push_back(classid);
	
    prepos = pos;
    preclassid = classid;
    while(A[pos][classid] >= 0)
    {                
	pos = A[prepos][preclassid];
	classid = C[prepos][preclassid];
	ybar.push_back(pos);//positions
	ybarlabel.push_back(classid);//class labels		
	prepos = pos;
	preclassid = classid;
    }
        
    marginloss = M[x.len][maxclassid];
    reverse(ybar.begin(), ybar.end());
    reverse(ybarlabel.begin(), ybarlabel.end());
}

/** Partial loss proposed by ref[2]
 *
 *  @param b1 [read] the newly selected starting position of a segment
 *  @param b2 [read] the newly selected starting position of the next segment
 *  @param y [read] the full position label i.e., all segment starting positions of a sequence
 *  @param ylabel [read] the full class label i.e., all segment class labels of a sequence
 *  @param classid [read ] the class id of segment [b1,b2)
 *  @param length [read ] the length of current sequence
 *  @return partial label loss value
 */
double CSMMMulticlassLoss::PartialDelta(const unsigned int b1, const unsigned int b2, 
	const vector<unsigned int> y, const vector<unsigned int> &ylabel, unsigned int classid, 
	unsigned int length)
{
    double segDelta = PartialSegDelta(b2,y);
    double labelDelta = PartialPointDelta(b1,b2,y,ylabel,classid);	
    double nearestSegDelta = PartialNearestSegDelta(b2,y,length);
    
    return (lossw[0]*segDelta + lossw[1]*labelDelta + lossw[2]*nearestSegDelta);
}


/** Partial segment-wise label loss proposed by ref[2]
 *  seeking for the nearest boundary y_i to ybar_i in boundary sequence y.
 *  delta = abs(y_i-ybar_i)
 * 
 *  @param ybar_i [read] the newly selected starting position of the next segment
 *  @param y [read] the full position label i.e., all segment starting positions of a sequence 
 *  @param length [read ] the length of current sequence
 *  @return partial label loss value
 */
double CSMMMulticlassLoss::PartialNearestSegDelta(const unsigned int ybar_i, 
	const vector<unsigned int> y, unsigned int length)
{	
    if(ybar_i>=length)
		return 0;
    if(binary_search(y.begin(), y.end(), ybar_i))
	return 0;
    unsigned int dist;
    unsigned int mindist = y[y.size()-1];
    for(unsigned int i = 0;i <y.size(); i++)
    {
	dist = (ybar_i>y[i]? (ybar_i-y[i]):(y[i]-ybar_i));		
	//cout << ybar_i<< "-" << y[i]<<" = " << dist << endl;
	if(mindist>dist)
	    mindist = dist;
    }
    //cout << "mindist = "<< mindist << endl;
    //getchar();
    
    return mindist;
}


/** Partial set-wise label loss proposed by ref[2]
 *
 *  @param ybar_i [read] the newly selected starting position of a segment
 *  @param y [read] the full label i.e., all segment starting positions of a sequence
 *
 *  @return partial label loss value
 */
double CSMMMulticlassLoss::PartialSegDelta(const unsigned int ybar_i, const vector<unsigned int> &y)
{
    if(binary_search(y.begin(), y.end(), ybar_i))
	return -1;
    else
	return 1;
}

/** Partial point-wise label loss proposed by ref[1]
 *
 *  @param b1 [read] the newly selected starting position of a segment
 *  @param b2 [read] the newly selected starting position of the next segment
 *  @param y [read] the full position label i.e., all segment starting positions of a sequence
 *  @param ylabel [read] the full class label i.e., all segment class labels of a sequence
 *  @param classid [read ] the class id of segment [b1,b2)
 *  @return partial label loss value
 */
double CSMMMulticlassLoss::PartialPointDelta(const unsigned int b1, const unsigned int b2, 
	vector<unsigned int> y, const vector<unsigned int> &ylabel, unsigned int classid)
{
    double delta = 0;
    vector<unsigned int>::iterator b1_up,b2_low;	
    b1_up = upper_bound(y.begin(),y.end(),b1);
    b2_low = lower_bound(y.begin(),y.end(),b2)-1;
    /** just need to compare region [b1,b1_up...b2_low,b2) including b1 and excluding b2
     */
    //printf("b1_up:%d,b2_low:%d",*b1_up,*b2_low);
    if((b1_up> b2_low) && (ylabel[b2_low-y.begin()]!=classid))
    {		
	delta+=b2-b1;
	return yeta*delta;
    }
    
    /* [b1,b1_up)*/
    int b1_up_pos = b1_up-y.begin();
    assert(b1_up_pos>0);
    if(ylabel[b1_up_pos-1]!=classid)
    {
	delta+=*b1_up-b1;
	//printf("*b1_up:%d,b1:%d\t",*b1_up,b1);
    }
    
    /* [b1_up...b2_low)*/
    int numOfSeg = b2_low-b1_up;
    int pos = b1_up_pos;
    for(int i=0;i<numOfSeg;i++)
    {
	if(ylabel[pos]!=classid)
	    delta+=y[pos+1]-y[pos];
	//printf("y[pos+1]:%d,y[pos]:%d\t",y[pos+1],y[pos]);
	pos++;
	
    }
    
    /* [b2_low...b2)*/
    unsigned int b2_low_pos = b2_low-y.begin();
    assert(b2_low_pos<y.size());
    if(ylabel[b2_low_pos]!=classid)
    {
	delta+=b2-*b2_low;
	//printf("b2:%d,*b2_low:%d\t",b2,*b2_low);
    }
    
    return yeta*delta;
}

double CSMMMulticlassLoss::Delta(const vector<unsigned int> &ybar, const vector<unsigned int> &y)
{
    unsigned int same = 0;
    unsigned int i=0, j=0;
    
    while(i<y.size() && j<ybar.size())
    {
	if(y[i] < ybar[j] && i<y.size()) 
	    i++;
	else if(y[i] > ybar[j] && j<ybar.size()) 
	    j++;
	else if(y[i] == ybar[j])
	{
	    same += 1;
	    i++;
	    j++;
	}
    }
    
    return (y.size()+ybar.size() - 2*same);
}

double CSMMMulticlassLoss::NearestSegDelta(const vector<unsigned int> &ybar, 
					   const vector<unsigned int> &y, unsigned int length)
{
    double sum = 0;
    unsigned int i;
    
    // from ybar to y
    for(i=0;i<ybar.size();i++)
    {
	sum += PartialNearestSegDelta(ybar[i],y,length);
    }        
    
    //from y to ybar
// 	for(i=0;i<y.size();i++)
// 	{
// 		sum += PartialNearestSegDelta(y[i],ybar,length);
// 	}        
    return sum;
}

double CSMMMulticlassLoss::AllDelta(const vector<unsigned int> &ybar, const vector<unsigned int> &y, 
	const vector<unsigned int> &ybarlabel, const vector<unsigned int> &ylabel, const int length)
{
	
    double segDelta,labelDelta,nearestSegDelta;
    
    vector <unsigned int> yss = Boundry2StatSequence(y,ylabel,length);
    vector <unsigned int> ybarss = Boundry2StatSequence(ybar,ybarlabel,length);
    
    segDelta = Delta(ybar,y);
    labelDelta = Labelloss(yss,ybarss);
    nearestSegDelta = NearestSegDelta(ybar,y,length);
    
    if(verbosity>=3)
	cout <<"seg_delta "<<segDelta << "\tlabel_delta "<<labelDelta<<"\tnearest_seg_delta "<<nearestSegDelta<<endl;
    
    return (lossw[0]*segDelta + lossw[1]*labelDelta + lossw[2]*nearestSegDelta);
}

        
 
/**  compute precision and recall
 *
 *   @param y [read] actual label
 *   @param ybar [read] predicted label
 *   @param prec [write] precision
 *   @param rec [write] recall
 */
void CSMMMulticlassLoss::PrecRec(const vector<unsigned int> &y, const vector<unsigned int> &ybar, double &prec, double &rec)
{
    double same = 0;
    prec = 0.0;
    rec = 0.0;
    
    
    //F_alpha = (1+alpha)*prec*rec/(alpha*prec+rec)
    unsigned int i=0, j=0;
    while(i<y.size() && j<ybar.size())
    {
	if(y[i] < ybar[j] && i<y.size()) 
	    i++;
	else if(y[i] > ybar[j] && j<ybar.size()) 
	    j++;
	else if(y[i] == ybar[j])
	{
	    same += 1;
	    i++;
	    j++;
	}
    }
    
    prec = same/ybar.size();
    rec = same/y.size();
}

/**  Evaluate the performance of the model on the examples in _data
 *
 *   @param model [read] model object containing a weight vector
 */

void CSMMMulticlassLoss::Evaluate(CModel* model)
{
    if(iscrossvalidation)
    {
       EvaluateAll(model,precision_all,recall_all,f1score_all,
                   allStd, allBar, allStd_pid,allPred_pid,acc_all_pid,acc_all_vote,acc_all);
       cout <<"Training Evaluation ----------------------"<< endl;
       EvaluateMark(model,SMM::TRAIN_DATA,precision_train,recall_train,f1score_train,
                    allStd_train,allBar_train,trainStd_pid,trainPred_pid,acc_train_pid,acc_train_vote,acc_train);
       cout <<"Testing Evaluation -----------------------"<< endl;
       EvaluateMark(model,SMM::TEST_DATA,precision_test,recall_test,f1score_test,
                    allStd_test,allBar_test,testStd_pid,testPred_pid,acc_test_pid,acc_test_vote,acc_test);
       cout <<"Validating Evaluation --------------------"<< endl;
       EvaluateMark(model,SMM::VALID_DATA,precision_valid,recall_valid,f1score_valid,
                    allStd_valid,allBar_valid,validStd_pid,validPred_pid,acc_validation_pid,acc_validation_vote,acc_validation);
       
    }
    else
       EvaluateAll(model,precision_all,recall_all,f1score_all,
                   allStd, allBar, allStd_pid,allPred_pid,acc_all_pid,acc_all_vote,acc_all);
}


/**  Evaluate the performance of the model on all the examples in _data
 *
 *   @param model [read] model object containing a weight vector
 */

void CSMMMulticlassLoss::EvaluateAll(CModel* model,double & weighted_prec,double & weighted_rec, 
                                     double& weighted_avgf1,vector <unsigned int> &yss,vector <unsigned int> &ybarss,
                                     vector<unsigned int> &pid,vector<unsigned int> & pid_bar,double & avgpid, 
                                     double & voted_acc, double & weighted_acc)
{
   // sanity check
   assert(model);
   assert(_data);
   assert(m > 0);
   
   yss.clear();
   ybarss.clear();
   pid.clear();
   pid_bar.clear();
   
   if(not _data->HasLabel())
   {
      throw CBMRMException("Data object does not contain labels","CSMMSegmentationLoss::Evaluate()");
   }
   
   // more defensive precaution
   assert(_data->HasLabel());
   assert(_data->dim() == model->GetDimOfW());
   
   TheMatrix &w = _model->GetW();
   const vector<CSeqMulticlassLabel::seqlabel_struct> &Y = _data->labels();
   const vector<CSeqMulticlassFeature::seqfeature_struct> &X = _data->features();
   
   double unweighted_avgf1 = 0;
   weighted_prec = 0;
   weighted_rec = 0;
   weighted_avgf1 = 0;
   double totalsegment = 0;
   weighted_acc = 0;
   unsigned int totalexamples = 0;
   voted_acc = 0;
   double correct_item_cnt = 0;
   int totalMarkeditems = 0;
   
   vector <int > cvmark = _data->Getcvmark();
   for(unsigned int i=0; i < m; i++)
   {
      pid.push_back(0);
      pid_bar.push_back(0);
      
      vector<unsigned int> ybar;
      vector<unsigned int> ybarlabel;
      double marginloss = 0;
      double marginloss_y = 0;
      double labelloss_y = 0;
      
      // find best label y' and return the score wrt to y'		
      if(is_single_action_persequence)
         find_best_label_grammer(X[i], w, ybar, ybarlabel, marginloss,0,_data->getNumOfClass());
      else
         find_best_label(X[i], w, ybar, ybarlabel, marginloss,0,_data->getNumOfClass());
      
      // accumulate loss                
      //double loss = marginloss + labelloss;                
      double prec = 0;
      double rec = 0;
      double f1 = 0;
      double acc = 0;
      
      PrecRec(Y[i].pos,ybar,prec,rec);  
      
      vector<unsigned int> yss_tmp;
      vector<unsigned int> ybarss_tmp;
      yss_tmp = Boundry2StatSequence(Y[i].pos,Y[i].type,X[i].len);
      ybarss_tmp = Boundry2StatSequence(ybar,ybarlabel,X[i].len);
      acc = Accuracy(yss_tmp,ybarss_tmp);
      
      for(unsigned int j=0; j < yss_tmp.size(); j++)
         yss.push_back(yss_tmp[j]);
      for(unsigned int j=0; j < ybarss_tmp.size(); j++)
         ybarss.push_back(ybarss_tmp[j]);
      
      f1 = (2*prec*rec)/(prec+rec);
      unweighted_avgf1 += f1;
      weighted_prec += Y[i].pos.size()*prec;
      weighted_rec += Y[i].pos.size()*rec;
      totalsegment += Y[i].pos.size();
      totalexamples += X[i].len;
      weighted_acc += X[i].len*acc;
      printf("ex: %2d   F1: %2.4f Acc:%2.4f\n",i,f1,acc);
      
      unsigned int voted_label = vote(ybarss_tmp);
      
      if(voted_label==yss_tmp[0])
         correct_item_cnt+=X[i].len;
      totalMarkeditems += X[i].len;
      
      
      if(verbosity>=1)
      {		
         printf("y:   ");
         for(unsigned int j=0;j<Y[i].pos.size();j++)
         {
            printf("%d(%d) ",Y[i].pos[j],Y[i].type[j]);
         }
         ComputeLoss(Y[i].pos,Y[i].type,Y[i].pos,Y[i].type,X[i],w,marginloss_y,labelloss_y,1);
         printf("marginloss:%f\n",marginloss_y);
         printf("ybar:");
         for(unsigned int j=0;j<ybar.size();j++)
         {
            printf("%d(%d) ",ybar[j],ybarlabel[j]);
         }
         printf("marginloss:%f\n",marginloss);
      }
   }
   avgpid = Accuracy(pid,pid_bar);
   weighted_prec /= totalsegment;
   weighted_rec /= totalsegment;
   weighted_avgf1 = 2*weighted_prec*weighted_rec/(weighted_prec+weighted_rec);
   unweighted_avgf1 /= m;
   weighted_acc /= totalexamples;
   voted_acc = correct_item_cnt / totalMarkeditems;
   
   cout << "1. unweighted macro-average F1: " << unweighted_avgf1 << endl;
   cout << "2. weighted macro-average F1  : " << weighted_avgf1 << endl;
   cout << "3. weighted accuracy          : " << weighted_acc << endl;        
   // micro or macro averaging f-measure
}

/**  Evaluate the performance of the model on all the test examples in _data
 *
 *   @param model [read] model object containing a weight vector
 */

void CSMMMulticlassLoss::EvaluateMark (CModel* model, int mark, double & weighted_prec,double & weighted_rec, 
				       double& weighted_avgf1,vector <unsigned int> &yss,vector <unsigned int> &ybarss,
				       vector<unsigned int> &pid,vector<unsigned int> & pid_bar,double & avgpid, 
				       double & voted_acc, double & weighted_acc)
{
    // sanity check
    assert(model);
    assert(_data);
    assert(m > 0);
    
    yss.clear();
    ybarss.clear();
    pid.clear();
    pid_bar.clear();
    
    if(mark!=SMM::TRAIN_DATA && mark!=SMM::TEST_DATA && mark!=SMM::VALID_DATA)
    {
	printf("mark:%d is invalid.\n",mark);
	exit(0);
    }
    
    if(not _data->HasLabel())
    {
	throw CBMRMException("Data object does not contain labels","CSMMSegmentationLoss::Evaluate()");
    }
    
    // more defensive precaution
    assert(_data->HasLabel());
    assert(_data->dim() == model->GetDimOfW());
    
    TheMatrix &w = _model->GetW();
    const vector<CSeqMulticlassLabel::seqlabel_struct> &Y = _data->labels();
    const vector<CSeqMulticlassFeature::seqfeature_struct> &X = _data->features();
    
    double unweighted_avgf1 = 0;
    weighted_prec = 0;
    weighted_rec = 0;
    weighted_avgf1 = 0;
    double totalsegment = 0;
    weighted_acc = 0;
    unsigned int totalexamples = 0;
    voted_acc = 0;
    double correct_item_cnt = 0;
    int totalMarkeditems = 0;
    
    vector <int > cvmark = _data->Getcvmark();
    
    for(unsigned int i=0; i < m; i++)
    {
	
	if(cvmark[i]!=mark)		
	    continue;
		
	pid.push_back(0);
	pid_bar.push_back(0);
	
	vector<unsigned int> ybar;
	vector<unsigned int> ybarlabel;
	double marginloss = 0;
	
	// find best label y' and return the score wrt to y'		
	if(is_single_action_persequence)
	    find_best_label_grammer(X[i], w, ybar, ybarlabel, marginloss,0,_data->getNumOfClass());
	else
	    find_best_label(X[i], w, ybar, ybarlabel, marginloss,0,_data->getNumOfClass());
	
	// accumulate loss                
	//double loss = marginloss + labelloss;                
	double prec = 0;
	double rec = 0;
	double f1 = 0;
	double acc = 0;
	
	PrecRec(Y[i].pos,ybar,prec,rec);  
	
	vector<unsigned int> yss_tmp;
	vector<unsigned int> ybarss_tmp;
	yss_tmp = Boundry2StatSequence(Y[i].pos,Y[i].type,X[i].len);
	ybarss_tmp = Boundry2StatSequence(ybar,ybarlabel,X[i].len);
	acc = Accuracy(yss_tmp,ybarss_tmp);
	
	for(unsigned int j=0; j < yss_tmp.size(); j++)
	    yss.push_back(yss_tmp[j]);
	for(unsigned int j=0; j < ybarss_tmp.size(); j++)
	    ybarss.push_back(ybarss_tmp[j]);
	
	f1 = (2*prec*rec)/(prec+rec);
	unweighted_avgf1 += f1;
	weighted_prec += Y[i].pos.size()*prec;
	weighted_rec += Y[i].pos.size()*rec;
	totalsegment += Y[i].pos.size();
	totalexamples += X[i].len;
	weighted_acc += X[i].len*acc;
	printf("ex: %2d   F1: %2.4f Acc:%2.4f\n",i,f1,acc);
	
	unsigned int voted_label = vote(ybarss_tmp);
	if(voted_label==yss_tmp[0])
	    correct_item_cnt+=X[i].len;
	totalMarkeditems += X[i].len;
	
	
	if(verbosity>=1)
	{		
	    printf("y:   ");
	    for(unsigned int j=0;j<Y[i].pos.size();j++)
	    {
		printf("%d(%d) ",Y[i].pos[j],Y[i].type[j]);
	    }
	    printf("\nybar:");
	    for(unsigned int j=0;j<ybar.size();j++)
	    {
		printf("%d(%d) ",ybar[j],ybarlabel[j]);
	    }
	    printf("marginloss:%f\n",marginloss);
	}
    }
    avgpid = Accuracy(pid,pid_bar);
    weighted_prec /= totalsegment;
    weighted_rec /= totalsegment;
    weighted_avgf1 = 2*weighted_prec*weighted_rec/(weighted_prec+weighted_rec);
    unweighted_avgf1 /= m;
    weighted_acc /= totalexamples;
    voted_acc = correct_item_cnt / totalMarkeditems;
    
    cout << "1. unweighted macro-average F1: " << unweighted_avgf1 << endl;
    cout << "2. weighted macro-average F1  : " << weighted_avgf1 << endl;
    cout << "3. weighted accuracy          : " << weighted_acc << endl;
    // micro or macro averaging f-measure
}


void CSMMMulticlassLoss::DecisionAndPrediction(CModel* model, int *ybar, double *f_ybar)
{}


vector <unsigned int> CSMMMulticlassLoss:: Boundry2StatSequence(const vector<unsigned int> segs,const vector<unsigned int> labels,unsigned int length)
{
    vector<unsigned int> statSequence;
    unsigned int i,j,end;
    for(i=0;i<segs.size();i++)
    {
	if(i==(segs.size()-1))
	    end = length-1;
	else
	    end = segs[i+1]-1;
	for(j=segs[i];j<=end;j++)
	{
	    statSequence.push_back(labels[i]);
	}
    }
    assert(statSequence.size()==length);
    return statSequence;
}

double CSMMMulticlassLoss:: Accuracy(const vector<unsigned int> target,const vector<unsigned int> predict)
{
    assert(target.size()>0&&predict.size()>0);
    unsigned int correctNum = 0;
    for(unsigned int i=0;i<target.size();i++)
    {
	if(target[i]==predict[i])		
	    correctNum++;
    }
    return (double)correctNum/(double)(target.size());
}

double CSMMMulticlassLoss:: Labelloss(const vector<unsigned int> target,const vector<unsigned int> predict)
{
    assert(target.size()>0&&predict.size()>0);
    unsigned int diffNum = 0;
    for(unsigned int i=0;i<target.size();i++)
    {
	if(target[i]!=predict[i])		
	    diffNum++;
    }
    return yeta*diffNum;
}

#endif
