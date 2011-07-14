#ifndef _COMMON_HPP_
#define _COMMON_HPP_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>



// #define's
// program mode
#define LEARNING_MODE      1    // do learning
#define PRED_AND_EVAL_MODE 2    // do prediction and evaluation (when labelset is supplied)
#define PREDICTION_MODE    3    // do prediction only (labelset is not supplied)

#ifdef PARALLEL_BMRM
// for use in parallel bmrm
#define ROOT_PROC_ID 0
#define MASTER(id) if((id) == ROOT_PROC_ID)
#define SLAVE(id)  if((id) != ROOT_PROC_ID)
#endif


// utility functions
bool IsBlankLine(std::string &line);

void trim(std::string& str);

void tokenize(const std::string& str, 
              std::vector<std::string>& tokens, 
              const std::string& delimiter = " ");


/** Functor for ascending indices sorting i.e., sort the indices based on the corresponding 
 *  values.  
 */
template<typename T>
struct indirect_less_than {
      T *array;
      indirect_less_than(T *value) : array(value) {}
      bool operator() (const int &a, const int &b) {return array[a] < array[b];}
};


/** Functor for descending indices sorting i.e., sort the indices based on the corresponding 
 *  values.  
 */
template<typename T>
struct indirect_greater_than {
      T *array;
      indirect_greater_than(T *value) : array(value) {}
      bool operator() (const int &a, const int &b) {return array[a] > array[b];}
};


/** Write array/vector data into file. Each element occupies a line.
 *
 *  @param fn [read] The file to store the data
 *  @param m [read] Number of elements in data
 *  @param v [read] The data
 */
template <class T>
bool WriteFile(const std::string& fn, const int& m, const T& v)
{
   using namespace std;
   ofstream fp(fn.c_str());
   
   if(not fp.good())
   {
      std::cout << "WriteFile(): Unable to open file (" + fn + ")" << std::endl;
      return false;
   }
   
   for(int i=0; i<m; i++)
      fp << v[i] << std::endl;
   
   fp.close();
   return true;
}



/** uniform pseudo random number generator (from NR)
 */
float ran1(long *idum);


/** generate pseudo random number from normal distribution (zero mean, unit variance) (from NR)
 */
float gasdev(long *idum);
#endif
