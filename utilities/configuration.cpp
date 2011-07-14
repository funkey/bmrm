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
 * Authors: Markus Weimer (markus.weimer@gmail.com)
 *          Choon Hui Teo (Choonhui.Teo@anu.edu.au)
 *
 * Created: (01/06/2007) 
 *
 * Last Updated: (28/10/2008)   
 */

#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "common.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

using namespace std;

void Configuration::SetDouble(string name, double value)
{
        IncreaseWriteCount(name);
        doubles[name] = value;  
}

void Configuration::SetDoubleVector(string name, vector<double> &value)
{
        IncreaseWriteCount(name);
        doublevectors[name] = value;  
}


double Configuration::GetDouble(string name) 
{  
        IncreaseReadCount(name);
        if(doubles.count(name)>0) return doubles[name];
        throw (CBMRMException("No double ("+name+") defined!",
                              "Configuration::GetDouble()"));  
}


vector<double>& Configuration::GetDoubleVector(string name) 
{  
        IncreaseReadCount(name);
        if(doublevectors.count(name)>0) return doublevectors[name];
        throw (CBMRMException("No double vector ("+name+") defined!",
                              "Configuration::GetDoubleVector()"));  
}


void Configuration::SetInt(string name, int value)
{  
        IncreaseWriteCount(name);
        ints[name] = value;  
}


void Configuration::SetIntVector(string name, vector<int> &value)
{  
        IncreaseWriteCount(name);
        intvectors[name] = value;  
}


int Configuration::GetInt(string name)
{  
        IncreaseReadCount(name);
        if(ints.count(name)>0) return ints[name];
        throw (CBMRMException("No int ("+name+") defined!",
                              "Configuration::GetInt()"));  
}


vector<int>& Configuration::GetIntVector(string name)
{  
        IncreaseReadCount(name);
        if(intvectors.count(name)>0) return intvectors[name];
        throw (CBMRMException("No int vector ("+name+") defined!",
                              "Configuration::GetIntVector()"));  
}



void Configuration::SetBool(string name, bool value)
{  
        IncreaseWriteCount(name);
        bools[name] = value;
}


void Configuration::SetBoolVector(string name, vector<bool> &value)
{  
        IncreaseWriteCount(name);
        boolvectors[name] = value;
}


bool Configuration::GetBool(string name)
{  
        IncreaseReadCount(name);
        if(bools.count(name)>0) return bools[name];
        throw (CBMRMException("No bool ("+name+") defined!",
                              "Configuration::GetBool()"));
  
}

vector<bool>& Configuration::GetBoolVector(string name)
{  
        IncreaseReadCount(name);
        if(boolvectors.count(name)>0) return boolvectors[name];
        throw (CBMRMException("No bool vector ("+name+") defined!",
                              "Configuration::GetBoolVector()"));
  
}



void Configuration::SetString(string name, string value)
{  
        strings[name] = value;
        IncreaseWriteCount(name);  
}


string Configuration::GetString(string name)
{
        IncreaseReadCount(name);
        if(strings.count(name)>0) return strings[name];
        throw (CBMRMException("No string ("+name+") defined!",
                              "Configuration::GetString()"));
}


/**   Check is a variable is set in the configuration file.
 *
 *    \param name [read] The name of the interested variable
 */
bool Configuration::IsSet(string name)
{  
        // check double
        if(doubles.count(name) > 0) return true;
        
        // check double vector
        if(doublevectors.count(name) > 0) return true;
        
        // check boolean
        if(bools.count(name) > 0) return true;
        
        // check boolean vector
        if(boolvectors.count(name) > 0) return true;
        
        // check integer
        if(ints.count(name) > 0) return true;               
        
        // check integer vector
        if(intvectors.count(name) > 0) return true;
        
        // check string
        if(strings.count(name) > 0) return true;          
        
        return false;
}


void Configuration::Dump(ostream& out)
{    
        out << "Usage statistics of configuration fields:" << endl;
        out << "type , name , readCount, writeCount" << endl;
  
        string s = " , ";
        
        // doubles
        for(d_itr iter = doubles.begin(); iter != doubles.end(); iter++)
        {
                const string& name = iter->first;
                out << "double , " << name << s << readCount[name] 
                    << s << writeCount[name] << endl;
        }
        
        // double vectors
        for(dvec_itr iter = doublevectors.begin(); iter != doublevectors.end(); iter++)
        {
                const string& name = iter->first;
                out << "doublevector , " << name << s << readCount[name] 
                    << s << writeCount[name] << endl;
        }
        
        // integers
        for(i_itr iter = ints.begin(); iter != ints.end(); iter++)
        {
                const string& name = iter->first;
                out << "int , " << name << s << readCount[name] 
                    << s << writeCount[name] << endl;
        }
        
        // integer vectors
        for(ivec_itr iter = intvectors.begin(); iter != intvectors.end(); iter++)
        {
                const string& name = iter->first;
                out << "intvector , " << name << s << readCount[name] 
                    << s << writeCount[name] << endl;
        }
        
        // booleans
        for(b_itr iter = bools.begin(); iter != bools.end(); iter++)
        {
                const string& name = iter->first;
                out << "bool , " << name << s << readCount[name] 
                    << s << writeCount[name] << endl;
        }
        
        // boolean vectors
        for(bvec_itr iter = boolvectors.begin(); iter != boolvectors.end(); iter++)
        {
                const string& name = iter->first;
                out << "boolvector , " << name << s << readCount[name] 
                    << s << writeCount[name] << endl;
        }
        
        // strings
        for(s_itr iter = strings.begin(); iter != strings.end(); iter++)
        {
                const string& name = iter->first;
                out << "string , " << name << s << readCount[name] 
                    << s << writeCount[name] << endl;
        }
}


void Configuration::IncreaseReadCount(const string& name)
{
        if(readCount.count(name) == 0)        
                readCount[name] = 1; 
        readCount[name]++;
}


void Configuration::IncreaseWriteCount(const string& name)
{  
        if(writeCount.count(name) == 0)
                writeCount[name] = 0;   
        writeCount[name]++;  
}


void Configuration::ReadFromFile(string filename)
{
        ifstream f(filename.c_str());
        if(not f.good())
        {
                throw CBMRMException("Unable to read configuration from " + filename, 
                                     "Configuration::ReadFromFile()");
        }
        
        while (! f.eof())
        {
                string line;
                getline (f,line);		
                trim(line);
                // Skip comments
                if(line.size() > 2 && line[0] == '/' && line[1] == '/') continue;
                // Skip blank lines
                if(line.size() == 0) continue;
    
                vector<string> tokens;
                tokenize(line, tokens, " ");
                if(tokens.size() < 3) continue;
    
                string name = tokens[1];
    
                if(tokens[0] == "double")
                {
                        double value = atof(tokens[2].c_str());
                        SetDouble(name, value);
                }
                else if(tokens[0] == "doublevector")
                {                        
                        vector<double> values;
                        for(size_t i=2; i<tokens.size(); i++)
                                values.push_back(atof(tokens[i].c_str()));
                        SetDoubleVector(name, values);
                }
                else if(tokens[0] == "int")
                {
                        int value = atoi(tokens[2].c_str());
                        SetInt(name, value);			
                }
                else if(tokens[0] == "intvector")
                {
                        vector<int> values;
                        for(size_t i=2; i<tokens.size(); i++)
                                values.push_back(atoi(tokens[i].c_str()));                        
                        SetIntVector(name, values);
                }
                else if(tokens[0] == "bool")
                {
                        trim(tokens[2]);
                        assert(tokens[2] == "true" || tokens[2] == "false"); 
                        string value = tokens[2];
                        trim(value);
                        if(value == "true") 
                                SetBool(name, true);
                        else
                        SetBool(name, false);
                }
                else if(tokens[0] == "boolvector")
                {
                        vector<bool> values;
                        for(size_t i=2; i<tokens.size(); i++)
                        {
                                trim(tokens[i]);
                                assert(tokens[i] == "true" || tokens[i] == "false"); 
                                string value = tokens[i];
                                if(value == "true") 
                                        values.push_back(true);
                                else
                                        values.push_back(false);
                        }        
                        SetBoolVector(name, values);
                }
                else if(tokens[0] == "string")
                {
                        string value = "";
                        for(size_t i = 2; i<tokens.size(); ++i) value += tokens[i] + " ";
                        trim(value);
                        SetString(name, value);
                }
                else
                {
                        throw CBMRMException("Cannot read line " + line + " from file " + filename, 
                                             "Configuration::ReadFromFile()");
                }
        }
        f.close();
}


void Configuration::Write(ostream& out)
{
        // doubles
        for(d_itr iter = doubles.begin(); iter != doubles.end(); iter++) 
        {
                out << "double " << iter->first << " " << iter->second << endl;
        }
        
        // double vectors
        for(dvec_itr iter = doublevectors.begin(); iter != doublevectors.end(); iter++) 
        {
                out << "doublevector " << iter->first << " ";
                vector<double> &values = iter->second;
                for(size_t i=0; i<values.size(); i++)
                        out << values[i] << " ";
                out << endl;                
        }
        
        // integers
        for(i_itr iter = ints.begin(); iter != ints.end(); iter++) 
        {
                out << "int " << iter->first << " " << iter->second << endl;
        }
        
        // integer vectors
        for(ivec_itr iter = intvectors.begin(); iter != intvectors.end(); iter++) 
        {
                out << "intvector " << iter->first << " ";
                vector<int> &values = iter->second;
                for(size_t i=0; i<values.size(); i++)
                        out << values[i] << " ";
                out << endl;                
        }
        
        // booleans
        for(b_itr iter = bools.begin(); iter != bools.end(); iter++) 
        {
                out << "bool " << iter->first << " " ;
                if (iter->second) 
                        out<< "true" << endl;
                else 
                        out<< "false" << endl;
        }
        
        // booleanr vectors
        for(bvec_itr iter = boolvectors.begin(); iter != boolvectors.end(); iter++) 
        {
                out << "boolvector " << iter->first << " ";
                vector<bool> &values = iter->second;
                for(size_t i=0; i<values.size(); i++)
                {
                        if(values[i])
                                out << "true ";
                        else
                                out << "false ";
                }                        
                out << endl;                
        }
        
        // strings
        for(s_itr iter = strings.begin(); iter != strings.end(); iter++) 
        {
                out << "string " << iter->first << " " << iter->second << endl;
        }        
}


void Configuration::WriteToFile(string fileName)
{  
        ofstream out(fileName.c_str());
        Write(out);
        out.close();
}


// Singleton implementation
Configuration* Configuration::instance = 0;

Configuration& Configuration::GetInstance(void)
{  
        if(instance == 0) 
                instance = new Configuration();
        return *instance;
}
