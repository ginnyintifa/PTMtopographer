//
//  input.hpp
//  
//
//  Created by Ginny Li on 19/10/15.
//
//

#ifndef input_feature_generation_hpp
#define input_feature_generation_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <tuple>
#include <cmath>
#include <cctype>

//http://www.boost.org/doc/libs/master/libs/algorithm/doc/html/algorithm/Searching.html


using namespace std;







struct Para_values
{
   
    int flank_length;
    //double stop_criterion;
    
    string fasta_name;
    //string peptide_name;
    
    string pred_sites;
    
    //double prob_threshold;
    
    string backinfo;
    
    string data_type;
    
    
};





vector<string> read_input(const string &fn)
{
    
    ifstream is(fn);
    check_init_istream(is);
    vector<string>myinput;
    
    string lines;
    
    for(;;)
    {
        getline(is,lines);
        
        if(!is) break;
        myinput.emplace_back(lines);
        
        
    }
    
    return myinput;
}


vector<string>get_parameterlines(vector<string>myinput)
{
    vector<string> parastring;
    for (size_t i=0; i<myinput.size(); i++)
    {
        if (myinput[i][0]=='>')
        {
            parastring.push_back(myinput[i]);
        }
        
        
    }
    
    return parastring;
    
}

//make changes later


Para_values  get_paravalues(vector<string>parastring)
{
    Para_values paravalues;
    
    
        for (size_t p=0; p<parastring[0].size(); p++)
        {
            if (parastring[0][p]=='=')
            {
                
                paravalues.flank_length=stoi(parastring[0].substr(p+1,(parastring[0].size()-p-1)));
                
                
            }
            
        }
    

    
    
    for (size_t p=0; p<parastring[1].size(); p++)
    {
        if (parastring[1][p]=='=')
        {
            
            paravalues.fasta_name=parastring[1].substr(p+1,(parastring[1].size()-p-1));
            
            
        }
        
    }

    

    for (size_t p=0; p<parastring[2].size(); p++)
    {
        if (parastring[2][p]=='=')
        {
            //just one site now
            
            paravalues.pred_sites=parastring[2].substr(p+1,(parastring[2].size()-p-1));
            
            
        }
        
    }

    
    
     for (size_t p=0; p<parastring[3].size(); p++)
    {
        if (parastring[3][p]=='=')
        {
           
            paravalues.backinfo=parastring[3].substr(p+1,(parastring[3].size()-p-1));
            
            
        }
        
    }
    
    
    
    for (size_t p=0; p<parastring[4].size(); p++)
    {
        if (parastring[4][p]=='=')
        {
            
            paravalues.data_type=parastring[4].substr(p+1,(parastring[4].size()-p-1));
            
            
        }
        
    }
    
    
    
    
        
    
    
    
    
       return paravalues;
}




#endif
