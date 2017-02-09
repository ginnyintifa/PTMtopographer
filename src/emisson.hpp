//
//  emisson.hpp
//  humansites
//
//  Created by Ginny Li on 15-8-25.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//

#ifndef humansites_emisson_hpp
#define humansites_emisson_hpp




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
#include <boost/algorithm/searching/boyer_moore_horspool.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

#include "math.hpp"





void assign_aaemissions(const std::vector<std::vector<std::pair<int, int>>>&states_windowrates,std::vector<std::vector<double>>&eachstate_20aa_emissons)
{
    using namespace std;
    
    for (size_t k=0; k<states_windowrates.size(); k++)
    {
        for (size_t p=0; p<states_windowrates[k].size(); p++)
        {
            eachstate_20aa_emissons[k].push_back(double(states_windowrates[k][p].first)/double(states_windowrates[k][p].second));
            
        }
        
        double sum=0;
        
        for (size_t o=0; o<eachstate_20aa_emissons[k].size(); o++)
        {
            sum=sum+eachstate_20aa_emissons[k][o];
        }
        
        
    }
    
   }


std::vector<std::vector<std::pair<double, double>>>get_emission(const std::vector<std::vector<double>>&hydro_2states)
{
	using namespace std;
	vector<vector<pair<double,double>>> emission_hydro{2};

    
    
	for (size_t h=0; h<2;h++)
    {
		emission_hydro[h].resize(hydro_2states[h].size());
        
		const double bandwidth=get_hvalue(hydro_2states[h]);
        
		vector<double>loggaussiankernel(hydro_2states[h].size());
        
  #pragma omp parallel for firstprivate(loggaussiankernel)
        
		for (size_t i=0;i<hydro_2states[h].size();i++)
        {
           // cout<<i<<endl;
            
			pair<double,double>hydro_prob;
			hydro_prob.first=hydro_2states[h][i];
			
			for (size_t p=0; p<hydro_2states[h].size(); p++)
            {
				loggaussiankernel[p]=loggaussian((hydro_2states[h][i]-hydro_2states[h][p])/bandwidth, 0, 1);
				
			}


			hydro_prob.second=logsuma(loggaussiankernel)-log(hydro_2states[h].size())-log(bandwidth);
			emission_hydro[h][i]=hydro_prob;
		}

	}
	return emission_hydro;
}




bool sortdouble(const double &adouble, const double &bdouble)
{
    return adouble<bdouble;
}


std::vector<std::vector<double>>choose_data(std::vector<std::vector<double>>&hydro_2states)
{
    using namespace std;
    
    
    vector<vector<double>>new_hydro2states{2};
    
    
    //I should choose data after sorting so that it is fair.
    
    //cout<<hydro_2states[0].size()<<"^^^^^^^^^^^^^^"<<endl;
    
    
    //cout<<hydro_2states[1].size()<<"^^^^^^^^^^^^"<<endl;
    
    //my aim is to uniformly choose 2000 data points
    
    //I think 2000 data points may not be enough let me try 5000 this time
    
    if (hydro_2states[0].size()>5000)
    {
        sort(hydro_2states[0].begin(),hydro_2states[0].end(),sortdouble);
        
        double themin=minimuma(hydro_2states[0]);
        double themax=maximuna(hydro_2states[0]);
       
        
        
        int sep=hydro_2states[0].size()/5000;//modified part
        

        
        for(size_t i=0;i<5000;i++)
        {
                new_hydro2states[0].push_back(hydro_2states[0][i*sep]);
            
                //cout<<hydro_2states[0][i*sep]<<" "<<i<<endl;
            
        }
        new_hydro2states[0].push_back(themin);
        new_hydro2states[0].push_back(themax);
        
        
       // cout<<"first one set**************"<<endl;
 
    }
    
    else  new_hydro2states[0]=hydro_2states[0];
    
    
    
    if (hydro_2states[1].size()>5000)
    {
        
        sort(hydro_2states[1].begin(),hydro_2states[1].end(),sortdouble);
        
        double themin1=minimuma(hydro_2states[1]);
       // cout<<themin1<<endl;
        double themax1=maximuna(hydro_2states[1]);
       // cout<<themax1<<endl;
        
   
        int sep1=hydro_2states[1].size()/5000;
        
        for(size_t i=0;i<5000;i++)
        {
           new_hydro2states[1].push_back(hydro_2states[1][i*sep1]);
           //cout<<hydro_2states[1][i*sep1]<<" "<<i<<endl;
            
        }
        new_hydro2states[1].push_back(themin1);
        new_hydro2states[1].push_back(themax1);
    }
    
    else  new_hydro2states[1]=hydro_2states[1];
    
    
        
    //cout<<"second one set*************"<<endl;

    return new_hydro2states;
    
    
}



bool sortvectorpair(const std::pair<double,std::vector<double>> &apair, const std::pair<double,std::vector<double>> &bpair)
{
    return apair.first<bpair.first;
}





std::vector<std::pair<double, std::vector<double>>>get_allemission(const std::vector<std::vector<double>>&hydro_forall,const std::vector<std::vector<std::vector<double>>>&hydro_allstates_emission)
{
    using namespace std;
    
    
    vector<pair<double,vector<double>>>hydro_allemisson;
    
    for (size_t i=0; i<hydro_forall.size(); i++)
    {
        for (size_t p=0; p<hydro_forall[i].size(); p++)
        {
            pair<double,vector<double>>each_block_emi;
            
            each_block_emi.first=hydro_forall[i][p];
            each_block_emi.second=hydro_allstates_emission[i][p];
            
            hydro_allemisson.push_back(each_block_emi);
        }
    }
    
    
    
    
    sort(hydro_allemisson.begin(),hydro_allemisson.end(),sortvectorpair);
    
    
    
    return hydro_allemisson;
    
    
}


#endif
