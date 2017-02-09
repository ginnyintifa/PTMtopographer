//
//  main.cpp
//  Ttest
//
//  Created by Ginny Li on 24/3/16.
//  Copyright © 2016 Ginny Li. All rights reserved.
//
#ifdef _MSC_VER
//#  pragma warning(disable: 4512) // assignment operator could not be generated.
//#  pragma warning(disable: 4510) // default constructor could not be generated.
//#  pragma warning(disable: 4610) // can never be instantiated - user defined constructor required.
#endif

#include <iostream>
#include <iomanip>
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


#include <boost/math/distributions/chi_squared.hpp>



using namespace std;

vector<vector<int>> cal_obs(vector<vector<int>>aa_compo_foreach,vector<int>ups,vector<int>downs)

{
    vector<vector<int>>chi_numbers_each{20};
    
    //cout<<"start to do the calculations..."<<endl;
    
    int ups_sum=0;
    int downs_sum=0;
    
    
    
    for (size_t u=0; u<ups.size(); u++)
    {
        for (size_t l=0; l<20; l++)
        {
            ups_sum=ups_sum+aa_compo_foreach[ups[u]][l];
        }
    }
    
    
    
    for (size_t d=0; d<downs.size(); d++)
    {
        for (size_t l=0; l<20; l++)
        {
            downs_sum=downs_sum+aa_compo_foreach[downs[d]][l];
        }
    }
    
  //  cout<<"got the summation..."<<endl;
    
    
        for (int l=0; l<20; l++)
        {
            
            for (size_t p=0; p<4; p++)
            {
                chi_numbers_each[l].push_back(0);

            }
            
           // cout<<"assigned the four conditions..."<<endl;
            
            for (size_t u=0; u<ups.size(); u++)
            {
                chi_numbers_each[l][0] =chi_numbers_each[l][0]+aa_compo_foreach[ups[u]][l];
            }
            
           // cout<<"get number 0 "<<chi_numbers_each[l][0]<<endl;
            
            chi_numbers_each[l][1]=ups_sum-chi_numbers_each[l][0];
            
            //cout<<"get number 1 "<<chi_numbers_each[l][1]<<endl;

            
            
            for (size_t d=0; d<downs.size(); d++)
            {
                chi_numbers_each[l][2] =chi_numbers_each[l][2]+aa_compo_foreach[downs[d]][l];

            }
            
            //cout<<"get number 2 "<<chi_numbers_each[l][2]<<endl;

            
            chi_numbers_each[l][3]=downs_sum-chi_numbers_each[l][2];
            
            //cout<<"get number 3 "<<chi_numbers_each[l][3]<<endl;


        }
    
    
   // cout<<"check an instance "<<chi_numbers_each[1][1]<<endl;
    
    
    return chi_numbers_each;
    
        
    
}











//in each protein, for each amino acid, find for group modified, how many T in this, and how many T in  non-modified group,


double cal_chi_stats(int modi_pos,int modi_neg, int nmodi_pos, int nmodi_neg )
{
    
  //  cout<<"get the chi statistics "<<endl;
    double chi_stats;
    
    int total=modi_pos+modi_neg+nmodi_pos+nmodi_neg;
    
    int modi_total=modi_pos+modi_neg;
    
    int nmodi_total=nmodi_pos+nmodi_neg;
    
    int pos_total=modi_pos+nmodi_pos;
    
    int neg_total=modi_neg+nmodi_neg;
    
    
    double modi_pos_exp= double(modi_total)*double(pos_total)/double(total);
    
    double modi_neg_exp= double(modi_total)*double(neg_total)/double(total);
    
    double nmodi_pos_exp= double(nmodi_total)*double(pos_total)/double(total);
    
    double nmodi_neg_exp= double(nmodi_total)*double(neg_total)/double(total);
    
    
    chi_stats=(modi_pos-modi_pos_exp)*(modi_pos-modi_pos_exp)/modi_pos_exp+(modi_neg-modi_neg_exp)*(modi_neg-modi_neg_exp)/modi_neg_exp+(nmodi_pos-nmodi_pos_exp)*(nmodi_pos-nmodi_pos_exp)/nmodi_pos_exp+(nmodi_neg-nmodi_neg_exp)*(nmodi_neg-nmodi_neg_exp)/nmodi_neg_exp;

//    cout<<"see what chi_stats is "<<chi_stats<<endl;

    return chi_stats;
    
    
}



double cal_pvalue(double chi_stats)
{
    boost::math::chi_squared mydist(1);
    
    //cout<<"calculate cdf..."<<endl;
    
    double pvalue=1-boost::math::cdf(mydist,chi_stats);
    
    //cout<<"my pvalue is... "<<pvalue<<endl;
    
    return pvalue;
}





