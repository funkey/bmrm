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
 * Authors: Markus Weimer       (markus.weimer@gmail.com)
 *          Choon Hui Teo       (ChoonHui.Teo@nicta.com.au)
 *
 * Created: 01/06/2007
 *
 * Last Updated: 10/01/2008
 */

#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

#include <map>
#include <string>
#include <cassert>
#include <vector>

#include "bmrmexception.hpp"

/**
 * Simple class to hold program configuration.
 *
 * The class is implemented as a singleton. To get a reference to the 
 * current instance, call Configuration::getConfiguration()
 */
class Configuration {
 private:
        void IncreaseReadCount(const std::string& name);
        void IncreaseWriteCount(const std::string& name);
  
        Configuration(void){}
        Configuration(const Configuration& other){assert(false);}

        typedef std::map<std::string,double> d_map;
        typedef std::map<std::string,double>::iterator d_itr;
        d_map doubles;
  
        typedef std::map<std::string,std::vector<double> > dvec_map;
        typedef std::map<std::string,std::vector<double> >::iterator dvec_itr;
        dvec_map doublevectors;
  
        typedef std::map<std::string,int> i_map;
        typedef std::map<std::string,int>::iterator i_itr;
        i_map ints;

        typedef std::map<std::string,std::vector<int> > ivec_map;
        typedef std::map<std::string,std::vector<int> >::iterator ivec_itr;
        ivec_map intvectors;
  
        typedef std::map<std::string,bool> b_map;
        typedef std::map<std::string,bool>::iterator b_itr;
        b_map bools;
  
        typedef std::map<std::string,std::vector<bool> > bvec_map;
        typedef std::map<std::string,std::vector<bool> >::iterator bvec_itr;
        bvec_map boolvectors;
  
        typedef std::map<std::string,std::string> s_map;
        typedef std::map<std::string,std::string>::iterator s_itr;
        s_map strings;
  
        typedef std::map<std::string, unsigned int> u_map;
        u_map readCount;
        u_map writeCount;
     
        static Configuration* instance;
  
  
public:
  
        /**
         * @return a reference to the current instance.
         */
        static Configuration& GetInstance(void);
  
        
        /**
         * Set a double value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetDouble(std::string name, double value);
  
        
        /**
         * Get a double value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	
        double GetDouble(std::string name);
  
        
        /**
         * Set a vector<double> value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetDoubleVector(std::string name, std::vector<double> &value);
  
        
        /**
         * Get a vector<double> value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	
        std::vector<double>& GetDoubleVector(std::string name);
        
        
        /**
         * Set a int value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetInt(std::string name, int value);
  
        
        /**
         * Get a int value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	
        int GetInt(std::string name);
        
        
        /**
         * Set a vector<int> value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetIntVector(std::string name, std::vector<int> &value);
  
        
        /**
         * Get a vector<int> value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	
        std::vector<int>& GetIntVector(std::string name);
        
        
        /**
         * Set a bool value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetBool(std::string name, bool value);
  
        
        /**
         * Get a bool value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	
        bool GetBool(std::string name);
  
        
        /**
         * Set a vector<bool> value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetBoolVector(std::string name, std::vector<bool> &value);
  
        
        /**
         * Get a vector<bool> value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	                
        std::vector<bool>& GetBoolVector(std::string name);
        
        
        /**
         * Set a string value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetString(std::string name, std::string value);
  
  
        /**
         * Get an string value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	
        std::string GetString(std::string name);
  
  
        /**
         * Set a vector<string> value in the configuration
         *
         * @param name the name of the configuration option.
         * @param value the value to be set.
         */
        void SetStringVector(std::string name, std::vector<std::string> &value);
  
  
        /**
         * Get an vector<string> value from the configuration
         *
         * @param name the name of the configuration option.
         * @throws CBMRMException if the given parameter cannot be found
         */	
        std::vector<std::string>& GetStringVector(std::string name);
        
        
        /** 
         * Is the given value set?
         * 
         * @param name 
         */
        bool IsSet(std::string name);
        
        
        /**
         * reads the configuration from file.
         *
         * The file format is very simple, each line holds one value. The lines have
         * the following format:
         *
         * [TYPE] [NAME] [VALUE]
         *
         * where [TYPE] can be: double, int, bool, or string
         *
         * Examples:
         *
         * int bmrm.maxIter 90
         * double bmrm.lambda 1.0
         * string logHeader This log was generated today. 
         * bool bias false
         * 
         * @param fileName the name of the file to read.
         */
        void ReadFromFile(std::string fileName);
        
  
        /**
         * Writes the config to the given stream.
         * 
         * This method uses the same format as described in readFromFile()
         *
         * @param out the stream where the configuration will be written to
         */
        void Write(std::ostream& out);
  
        
        /**
         * Writes the configuration to a file.
         *
         * This is a convenience wrapper around write()
         *
         * @param fileName the name of the file where the config will be stored.
         * @see write()
         */
        void WriteToFile(std::string fileName);
  
  
        /**
         * Dumps the usage statistics for all known fields to the given 
         * ostream. 
         *
         * Very useful for debugging.
         * 
         * @param out the stream to write to.
         */
        void Dump(std::ostream& out);
};

#endif /* _CONFIGURATION_HPP_ */
