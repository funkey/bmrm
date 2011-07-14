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
 *
 * Created: (14/02/2008)
 *
 * Last Updated: (14/02/2008)
 */

#ifndef _LINESEARCH_HPP_
#define _LINESEARCH_HPP_


#include <algorithm>
#include <stack>
#include <iterator>
#include <queue>
#include <vector>
#include "sml.hpp"
#include "model.hpp"
#include "loss.hpp"
#include "data.hpp"
#include "timer.hpp"
//#include "linesearch_utilities.hpp"
#define MYZERO 1e-10

/**
 * The basic class for line search
 */
class CLinesearch {
      
   protected:
      /**
       * Step size
       */
      double eta;
      
      /**
       * Pointer to the model object
       */
      CModel* _model;
      
      /**
       * Number of training instances
       */
      int numElements;
      
      /**
       * Number of classes
       */
      int numClass;
      
      /**
       * Scaling factor for the sum of the losses,
       * by default it is set to be 1/numElements
       */
      double C;
      
      /**
       * The sum of all losses divided by the number of elements
       */
      double error;
      
   public:
      /**
       * Class constructor.
       *
       * @param model [r/w] pointer to the model object
       * @param data [r/w] pointer to the loss object
       */
      CLinesearch(CModel* model, CLoss* loss, CData *data);
      
      /**
       * Destructor
       */
      virtual ~CLinesearch() {}
      
      
      /**
       * Return the averaged sum of losses.
       */
      virtual double GetLoss() 
      {
         return error;
      }
};



/*
  LINESEARCH UTILITIES
 */

/** for sorting indices */
class ScalarWithIndex {
   public:
      double value;
      int index;
      ScalarWithIndex(double val, int ind) :
         value(val), index(ind) {
      }
      
      ScalarWithIndex() :
         value(0.0), index(-1) {
      }
      
      ~ScalarWithIndex() {
      }
      
      inline int operator<(const ScalarWithIndex & other) const {
         return (value < other.value);
      }
      
      inline void swap(ScalarWithIndex & other) {
         double myvalue = value;
         int myindex = index;
         value = other.value;
         index = other.index;
         other.value = myvalue;
         other.index = myindex;
      }
};

/**
 * Class that defines a pair of labels
 */
class LabelPair {
   public:
      int trueLabel;
      int label;
      
      LabelPair(int z, int z_prime) {
         trueLabel = z;
         label = z_prime;
      }
      
      LabelPair() :
         trueLabel(-1), label(-1) {
      }
      
      ~LabelPair() {
      }
};

/**
 * Class that defines a linear line
 */
class Line {
   public:
      double value;
      LabelPair labels;
      double slope;
      double bias;
      
      Line(double val, LabelPair z_zprime, double s, double b) :
         value(val), labels(z_zprime), slope(s), bias(b) {
      }
      
      
      Line() :
         value(0.0), labels(LabelPair(-1, -1)), slope(0), bias(0) {
      }
      
      
      ~Line() {
      }
      
      
      inline int operator<(const Line & other) const {
         return (value < other.value);
      }
      
      inline void swap(Line & other) {
         double myValue = value;
         LabelPair mylabels = labels;
         double mySlope = slope;
         double myBias = bias;
         
         value = other.value;
         labels = other.labels;
         slope = other.slope;
         bias = other.bias;
         
         other.value = myValue;
         other.labels = mylabels;
         other.slope = mySlope;
         other.bias = myBias;
      }
};

/**
 * A structure that defines a linear line with 'step' information
 */
struct lineWithLabel {
      double step;
      double slope;
      double bias;
      std::vector<LabelPair> labels;
};


/**
 * Class that defines a pair of labels that is associated with an index
 */
class LabelWithIndex {
   public:
      std::vector<LabelPair> labels;
      int index;
      
      LabelWithIndex(const std::vector<LabelPair> &l, const int& ind) {
         labels = l;
         index = ind;
      }
      
      virtual ~LabelWithIndex() {
      }
};

/**
 * Class that defines a linear line that is associated with some labels and an index
 */
class lineWithLabelAndIndex {
   public:
      std::deque<lineWithLabel>::iterator line;
      int index;
      
      lineWithLabelAndIndex(const std::deque<lineWithLabel>::iterator &l,
                            const int& ind) {
         line = l;
         index = ind;
      }
      
      virtual ~lineWithLabelAndIndex() {
      }
      
      /**
       * Define the greater operation
       * @param a [read] the object that we compare with
       */
      bool operator>(const lineWithLabelAndIndex &a) const {
         return line->step > (a.line)->step;
      }
};


/**
 * Class that uses a priority queue to pick the line that has the
 * minimal 'step' value. All lines are stored in a vector
 * of double-ended queues; and we assume line saved in each double-ended
 * queue is already sorted according to their 'step' values.
 */
class CVector_Operations {
   protected:
      /**
       * A priority queue, having the line with the minimal 'step' value on the top
       */
      std::priority_queue<lineWithLabelAndIndex, std::vector<lineWithLabelAndIndex>, std::greater<lineWithLabelAndIndex> > minQ;

   public:
      /**
       * Class constructor
       */
      CVector_Operations(std::vector<std::deque<lineWithLabel> >& vec_deque, std::vector<std::deque<
                         lineWithLabel>::iterator> &vec_iterator) {
         int ind = 0;
         std::vector<std::deque<lineWithLabel> >::iterator vec_deque_it =
            vec_deque.begin();
         for (std::vector<std::deque<lineWithLabel>::iterator>::iterator it =
                 vec_iterator.begin(); it < vec_iterator.end(); it++, vec_deque_it++) {
			if (*it >= vec_deque_it->begin()) {
               lineWithLabelAndIndex lineWithInd(*it, ind);
               minQ.push(lineWithInd);
               (*it)--;
			}
			ind++;
         }
      }

      /**
       * Default class destructor
       */
      virtual ~CVector_Operations() {
      }
      
      /**
       * Return the priority queue
       */
      virtual std::priority_queue<lineWithLabelAndIndex,std::vector<lineWithLabelAndIndex> , std::greater<lineWithLabelAndIndex> > getQ() {
         return minQ;
      }
      
      /**
       * Return a vector of lines that has the minimal 'step' value
       * @param vec_deque [read/write] linear lines that are stored in a vector of double-ended queue
       * @param vec_iterator [read/write] iterator of lines that we will consider in the next iterate
       */
      virtual std::vector<lineWithLabelAndIndex> Min(
         std::vector<std::deque<lineWithLabel> >& vec_deque, std::vector<std::deque<
         lineWithLabel>::iterator> &vec_iterator) {
         std::vector<lineWithLabelAndIndex> minimal;
         if (!minQ.empty()) {
			minimal = std::vector<lineWithLabelAndIndex> (1, minQ.top());
			double minVal = (minQ.top().line)->step;
			minQ.pop();
            
			//pop up every line with the same step as the minimal step
			while ((!minQ.empty()) && (minQ.top().line)->step == minVal) {
               minimal.push_back(minQ.top());
               minQ.pop();
			}
            
			//push in new lines
			for (std::vector<lineWithLabelAndIndex>::iterator it = minimal.begin(); it
					< minimal.end(); it++) {
               std::deque<lineWithLabel>::iterator itt = vec_iterator.at(it->index);
               if (itt >= (vec_deque.at(it->index)).begin()) {
                  lineWithLabelAndIndex newLine(itt, it->index);
                  (vec_iterator.at(it->index))--;
                  minQ.push(newLine);
               }
			}
         }
         return minimal;
      }
};

#endif
