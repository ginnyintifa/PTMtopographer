//
//  input.hpp
//  
//
//  Created by Ginny Li on 19/10/15.
//
//

#ifndef input_additional_annotation_hpp
#define input_additional_annotation_hpp

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
    string protid;
    string domain_start;
    string domain_end;
    string domain_name;
    string domain_type;
    
    string can_sites_position;
    string prots_id;
    
    string prot_cytoskeleton;
    string prot_cytosol;
    string prot_endoplasmic;
    string prot_endosome;
    string prot_extracellular;
    string prot_golgi;
    string prot_lysosome;
    string prot_membrane;
    string prot_mitochondrion;
    string prot_nucleus;
    string prot_peroxisome;
    
    string rf_score;
    
    

    string head_annotation;
    
    string global_dfdr;
    
    string prot_specific_dfdr;
    
    string prot_head_annotation;
    
    
    
    double first_fdr_score;
    double second_fdr_score;
    double third_fdr_score;
   
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
                
                paravalues.protid=parastring[0].substr(p+1,(parastring[0].size()-p-1));
                
                
            }
            
        }
    

    
    
    for (size_t p=0; p<parastring[1].size(); p++)
    {
        if (parastring[1][p]=='=')
        {
            
            paravalues.domain_start=parastring[1].substr(p+1,(parastring[1].size()-p-1));
            
            
        }
        
    }

    

    for (size_t p=0; p<parastring[2].size(); p++)
    {
        if (parastring[2][p]=='=')
        {
            //just one site now
            
            paravalues.domain_end=parastring[2].substr(p+1,(parastring[2].size()-p-1));
            
            
        }
        
    }

    
    
     for (size_t p=0; p<parastring[3].size(); p++)
    {
        if (parastring[3][p]=='=')
        {
           
            paravalues.domain_name=parastring[3].substr(p+1,(parastring[3].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[4].size(); p++)
    {
        if (parastring[4][p]=='=')
        {
            
            paravalues.domain_type=parastring[4].substr(p+1,(parastring[4].size()-p-1));
            
            
        }
        
    }

    
    
    for (size_t p=0; p<parastring[5].size(); p++)
    {
        if (parastring[5][p]=='=')
        {
            
            paravalues.can_sites_position=parastring[5].substr(p+1,(parastring[5].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[6].size(); p++)
    {
        if (parastring[6][p]=='=')
        {
            
            paravalues.prots_id=parastring[6].substr(p+1,(parastring[6].size()-p-1));
            
            
        }
        
    }
    


    for (size_t p=0; p<parastring[7].size(); p++)
    {
        if (parastring[7][p]=='=')
        {
            
            paravalues.prot_cytoskeleton=parastring[7].substr(p+1,(parastring[7].size()-p-1));
            
            
        }
        
    }
    for (size_t p=0; p<parastring[8].size(); p++)
    {
        if (parastring[8][p]=='=')
        {
            
            paravalues.prot_cytosol=parastring[8].substr(p+1,(parastring[8].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[9].size(); p++)
    {
        if (parastring[9][p]=='=')
        {
            
            paravalues.prot_endoplasmic=parastring[9].substr(p+1,(parastring[9].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[10].size(); p++)
    {
        if (parastring[10][p]=='=')
        {
            
            paravalues.prot_endosome=parastring[10].substr(p+1,(parastring[10].size()-p-1));
            
            
        }
        
    }
    for (size_t p=0; p<parastring[11].size(); p++)
    {
        if (parastring[11][p]=='=')
        {
            
            paravalues.prot_extracellular=parastring[11].substr(p+1,(parastring[11].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[12].size(); p++)
    {
        if (parastring[12][p]=='=')
        {
            
            paravalues.prot_golgi=parastring[12].substr(p+1,(parastring[12].size()-p-1));
            
            
        }
        
    }
    for (size_t p=0; p<parastring[13].size(); p++)
    {
        if (parastring[13][p]=='=')
        {
            
            paravalues.prot_lysosome=parastring[13].substr(p+1,(parastring[13].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[14].size(); p++)
    {
        if (parastring[14][p]=='=')
        {
            
            paravalues.prot_membrane=parastring[14].substr(p+1,(parastring[14].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[15].size(); p++)
    {
        if (parastring[15][p]=='=')
        {
            
            paravalues.prot_mitochondrion=parastring[15].substr(p+1,(parastring[15].size()-p-1));
            
            
        }
        
    }
    for (size_t p=0; p<parastring[16].size(); p++)
    {
        if (parastring[16][p]=='=')
        {
            
            paravalues.prot_nucleus=parastring[16].substr(p+1,(parastring[16].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[17].size(); p++)
    {
        if (parastring[17][p]=='=')
        {
            
            paravalues.prot_peroxisome=parastring[17].substr(p+1,(parastring[17].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[18].size(); p++)
    {
        if (parastring[18][p]=='=')
        {
            
            paravalues.rf_score=parastring[18].substr(p+1,(parastring[18].size()-p-1));
            
            
        }
        
    }
     
    for (size_t p=0; p<parastring[19].size(); p++)
    {
        if (parastring[19][p]=='=')
        {
            
            paravalues.head_annotation=parastring[19].substr(p+1,(parastring[19].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[20].size(); p++)
    {
        if (parastring[20][p]=='=')
        {
            
            paravalues.global_dfdr=parastring[20].substr(p+1,(parastring[20].size()-p-1));
            
            
        }
        
    }
    
    for (size_t p=0; p<parastring[21].size(); p++)
    {
        if (parastring[21][p]=='=')
        {
            
            paravalues.prot_specific_dfdr=parastring[21].substr(p+1,(parastring[21].size()-p-1));
            
            
        }
        
    }
    
    
    
    for (size_t p=0; p<parastring[22].size(); p++)
    {
        if (parastring[22][p]=='=')
        {
            
            paravalues.prot_head_annotation=parastring[22].substr(p+1,(parastring[22].size()-p-1));
            
            
        }
        
    }
    
    
    for (size_t p=0; p<parastring[23].size(); p++)
    {
        if (parastring[23][p]=='=')
        {
            
            paravalues.first_fdr_score=stod(parastring[23].substr(p+1,(parastring[23].size()-p-1)));
            
            
        }
        
    }
    
    
    
    for (size_t p=0; p<parastring[24].size(); p++)
    {
        if (parastring[24][p]=='=')
        {
            
            paravalues.second_fdr_score=stod(parastring[24].substr(p+1,(parastring[24].size()-p-1)));
            
            
        }
        
    }
    
    
    
    for (size_t p=0; p<parastring[25].size(); p++)
    {
        if (parastring[25][p]=='=')
        {
            
            paravalues.third_fdr_score=stod(parastring[25].substr(p+1,(parastring[25].size()-p-1)));
            
            
        }
        
    }
    
    
    

    
    
    
    
    
       return paravalues;
}




#endif
