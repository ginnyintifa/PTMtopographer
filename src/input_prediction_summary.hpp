//
//  input.hpp
//  
//
//  Created by Ginny Li on 19/10/15.
//
//

#ifndef input_prediction_summary_hpp
#define input_prediction_summary_hpp

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
// #include "boyer_moore_horspool.hpp"
//#include <boost/algorithm/searching/boyer_moore_horspool.hpp>
//#include <boost/foreach.hpp>
//#include <boost/functional/hash.hpp>



using namespace std;

/*
struct Para_values
{
    int block_length;
    int block_overlaps;
    double stop_criterion;
};
*/

//future version



struct Para_values
{
    string pred_probs;
    
    string decoy_probs;
    
    string my_prots;
    
    string can_prots;
    
    string decoy_prots;
    
    //string pred_states;
    
    int can_col;
    int decoy_col;
    
    string decoy_tag;
    
  };



/*
template <class CharT, class Traits>
void check_init_istream(::basic_istream<CharT, Traits>& is)
{
    is.exceptions(::ios_base::failbit);
    is.exceptions(::ios_base::badbit);
}
*/


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
            
            paravalues.pred_probs=parastring[0].substr(p+1,(parastring[0].size()-p-1));
            
            
        }
        
    }

    
    for (size_t p=0; p<parastring[1].size(); p++)
    {
        if (parastring[1][p]=='=')
        {
            
            paravalues.decoy_probs=parastring[1].substr(p+1,(parastring[1].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[2].size(); p++)
    {
        if (parastring[2][p]=='=')
        {
            
            paravalues.my_prots=parastring[2].substr(p+1,(parastring[2].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[3].size(); p++)
    {
        if (parastring[3][p]=='=')
        {
            
            paravalues.can_prots=parastring[3].substr(p+1,(parastring[3].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[4].size(); p++)
    {
        if (parastring[4][p]=='=')
        {
            
            paravalues.decoy_prots=parastring[4].substr(p+1,(parastring[4].size()-p-1));
            
            
        }
        
    }
    
    /*
    
    for (size_t p=0; p<parastring[5].size(); p++)
    {
        if (parastring[5][p]=='=')
        {
            
            paravalues.pred_states=parastring[5].substr(p+1,(parastring[5].size()-p-1));
            
            
        }
        
    }
    
    
    */
    
    
    for (size_t p=0; p<parastring[5].size(); p++)
    {
        if (parastring[5][p]=='=')
        {
            
            paravalues.can_col=stoi(parastring[5].substr(p+1,(parastring[5].size()-p-1)));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[6].size(); p++)
    {
        if (parastring[6][p]=='=')
        {
            
            paravalues.decoy_col=stoi(parastring[6].substr(p+1,(parastring[6].size()-p-1)));
            
            
        }
        
    }

    
    
    
    for (size_t p=0; p<parastring[7].size(); p++)
    {
        if (parastring[7][p]=='=')
        {
            
            paravalues.decoy_tag=parastring[7].substr(p+1,(parastring[7].size()-p-1));
            
            
        }
        
    }
    

    
    
       return paravalues;
}


















#endif
