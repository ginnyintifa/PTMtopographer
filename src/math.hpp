//
//  math.hpp
//  humansites
//
//  Created by Ginny Li on 15-8-30.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//

#ifndef math_hpp
#define math_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
//#include <unordered_set>
#include <algorithm>
#include <limits>
//#include <tuple>
#include <cstdio>
#include <cctype>
//http://www.boost.org/doc/libs/master/libs/algorithm/doc/html/algorithm/Searching.html
//#include <boost/algorithm/searching/boyer_moore_horspool.hpp>
//#include <boost/functional/hash.hpp>
#include <cmath>

#include <sstream>


using namespace std;


double large_num=4;

//gaussian density

double gaussianpdf(double ratio, double mean, double variance)
{
    double prob;
    prob=(1/(sqrt(2*M_PI*variance)))*exp( -0.5*pow((ratio-mean),2.0)/variance);
    
    return prob;
    
}

//log gaussian density

double loggaussian(double ratio, double mean, double variance)
{
    //double logprob;
    
    return -(ratio-mean)*(ratio-mean)/(2.0*variance)-log(sqrt(2.0*M_PI*variance));
    
    //return logprob;
}


//get maximum

//now I should have the function for 4

//for these functions below I should unify them in the future
double maximum(double va, double vb)
{
    double max;
    max=va;
    
    if (vb>max)
    {
        max=vb;
    }
    
        return max;
}



//arg maximum
int maxarg(double va, double vb)
{
    int num;
    double max;
    
    max=va;
    num=0;
    
    if (vb>max)
    {
        max=vb;
        num=1;
    }
    
   
    
        return num;
}


//mean

double getmean(vector<double>ratios)
{
    double sumratio=0;
    for (size_t i=0;i<ratios.size();i++)
    {
        sumratio=sumratio+ratios[i];
    }
    
    return double(sumratio)/double(ratios.size());
}


//sum of integers
int sumofintegers(vector<int>signs)
{
    int sum=0;
    for(size_t i=0;i<signs.size();i++)
    {
        sum=sum+signs[i];
    }
    return sum;
    
}


//variance

double getvariance(vector<double>ratios,double mean)
{
    double sumsq=0;
    for (size_t i=0;i<ratios.size();i++)
    {
        sumsq=sumsq+(ratios[i]-mean)*(ratios[i]-mean);
    }
    return double(sumsq)/double(ratios.size());
}



//log scale summation

//want to get the log of summatin of A, B but you only know a, b, where a=logA, b=logB



double logsum(double firstv,double secondv)
{
    double maxv=maximum(firstv, secondv);
    
    double nfirstv=firstv-maxv;
    double nsecondv=secondv-maxv;
   // double nthirdv=thirdv-maxv;
    //double nfourdv=fourdv-maxv;
    
    return maxv+log(exp(nfirstv)+exp(nsecondv));
    
}







//maximum of 9 numbers

double max4(double logt[4])
{
    double maxt=logt[0];
    for (int i=0; i<4;i++)
    {
        if (logt[i]>maxt)
            maxt=logt[i];
        
    }
    
    return maxt;
}


//log scale summation of an array of 9 numbers

double logsum4(double logt[4])
{
    double exp4sum=0;
    for (int i=0; i<4;i++)
    {
        // logt[i]-max9(logt);
        
        exp4sum=exp4sum+exp(logt[i]-max4(logt));
    }
    
    return (max4(logt)+log(exp4sum));
    
    
}




double maximuna(vector<double> nums)
{
    int c=nums.size();
    double maxnum=nums[c-1];
    for (int i=0; i<nums.size();i++)
    {
        if (nums[i]>maxnum)
            maxnum=nums[i];
    }
    return maxnum;
    
    
}



double minimuma(vector<double> nums)
{
    int c=nums.size();
    double minnum=nums[c-1];
    for (int i=0; i<nums.size();i++)
    {
        if (nums[i]<minnum)
            minnum=nums[i];
    }
    return minnum;
    
    
}

int get_argminimuma(vector<double> nums)
{
    int c=nums.size();
    double minnum=nums[c-1];
    int argmin=c-1;
    for (int i=0; i<nums.size();i++)
    {
        if (nums[i]<minnum)
        {
            minnum=nums[i];
            argmin=i;
        }
    }
    return argmin;
    
}


int get_argmaximuma(vector<double> nums)
{
    int c=nums.size();
    double maxnum=nums[c-1];
    int argmax=c-1;
    for (int i=0; i<nums.size();i++)
    {
        if (nums[i]>maxnum)
        {
            maxnum=nums[i];
            argmax=i;
        }
    }
    return argmax;
    
}


//this function can not calculate too big vectors
/*double logsuma(vector<double> logt)
{
    double expsum=0;
    double max=maximuna(logt);
    for(size_t i=0;i<logt.size();i++)
    {
        expsum=expsum+exp(logt[i]-max);
    }
    
    return max+log(expsum);
    
}
*/




//simple summation

double logsumsimple(double firstn,double secondn, double thirdn)
{
    
    return log(exp(firstn)+exp(secondn)+exp(thirdn));
}


//read in...






//read in multiple sequences of observations






vector<vector<double>> read_tsv_manylines_double(const string fn)
{
    
    ifstream myfile (fn);
    
    vector<double> everyline;
    vector<vector<double>>mydoubles;
    
    string line;
    
    while(getline(myfile,line))
    {
        stringstream stream(line);
        double ratio;
        while(stream>>ratio)
        {
            everyline.push_back(ratio);
        }
        mydoubles.push_back(everyline);
        everyline.clear();
    }
    
    
    return mydoubles;
}




//a function to calculate probability from multinomial distribution

//calculate on log scale


//lets see what is wrong here


double get_logmultinomial(vector<int>myaa_obs,vector<double>eachstate_20aa_emi)// for a specific state one result
{
    //double theresult;
    
    
    size_t mylen=20;//myaa_obs.size();//suppose to be the vector size of each block 20//yes I ignored the last 4 possibilities
    
   // cout<<"test the length "<<mylen<<endl;
    
    double logmultinomialc=0;
    double logmultinomialm=0;
    
    double logmultinomial;
    
    double blocksize=0;
    
    for (size_t i=0; i<mylen; i++)
    {
        blocksize+=myaa_obs[i];
        
        
        for (int h=1; h<myaa_obs[i]+1; h++)
        {
            logmultinomialc+=log(h);
        }
        
        logmultinomialm+=log(eachstate_20aa_emi[i])*myaa_obs[i];
        
    }
    
    double mylenf=0;
    
    for (int u=1;u<blocksize+1;u++)
    {
        mylenf+=log(u);
    }
        
    
    
    logmultinomialc=mylenf-logmultinomialc;
    
    
    
    logmultinomial=logmultinomialc+logmultinomialm;
    
    
    
    return logmultinomial;
    
    
    
}


//I think I should first arrange my emissions in order accenting?
//in the first phase I don't need to order them actually

//now a function to give the emission probability once I have a observation of known state and value//and they should be in log scale


//so How do i interpolate? this is a question,I want to use the closest to approximate,yes this is the first phase

double get_nonpara(double obsvalue,int astate,const vector<vector<pair<double,double>>>& aproperty_emission)
{
    
    
    vector<double>distance;
    
    for (size_t i=0;i<aproperty_emission[astate].size();i++)
    {
        distance.push_back(abs(obsvalue-aproperty_emission[astate][i].first));
        
    }
    
    int argclosest=get_argminimuma(distance);
    
    double thedensity=aproperty_emission[astate][argclosest].second;
    
    return thedensity;
    
}



//the above is the simplest interpolate so what I will do is to implement a linear interpolation


double get_lnonpara(double obsvalue,int astate,const vector<vector<pair<double,double>>>& aproperty_emission)
{
    
    
    double thedensity=0;//new modi
    
    pair<vector<double>,vector<int>>distancepos;
    pair<vector<double>,vector<int>>distanceneg;
    
    
    
    for (size_t i=0;i<aproperty_emission[astate].size();i++)
    {
        pair<double, int> adis;
        
        adis.first=obsvalue-aproperty_emission[astate][i].first;
        adis.second=i;
        
        if (adis.first>0)
        {
            
            
            distancepos.first.push_back(adis.first);
            distancepos.second.push_back(adis.second);
            
        }
        else
        {
            distanceneg.first.push_back(adis.first);
            distanceneg.second.push_back(adis.second);
        }
        
    }
    
    //i think there is something wrong here
    //yes, need to change
    
    
    
    if(distancepos.first.size()>0&&distanceneg.first.size()>0)
    {
    
    int argclosestpos=get_argminimuma(distancepos.first);
    int argclosestneg=get_argmaximuma(distanceneg.first);
        
        
        //may have coredumped in this area
    
    double xpos=aproperty_emission[astate][distancepos.second[argclosestpos]].first;
    double xneg=aproperty_emission[astate][distanceneg.second[argclosestneg]].first;

    
    double denpos=aproperty_emission[astate][distancepos.second[argclosestpos]].second;
    double denneg=aproperty_emission[astate][distanceneg.second[argclosestneg]].second;
    
    
    
    
     thedensity=(denpos-denneg)*(obsvalue-xneg)/(xpos-xneg)+denneg;
    }
    
    else if (distancepos.first.size()>0)
    {
        int argpos=get_argminimuma(distancepos.first);
        thedensity=aproperty_emission[astate][distancepos.second[argpos]].second;

    }
    
    else if (distanceneg.first.size()>0)
    {
        int argneg=get_argmaximuma(distanceneg.first);
        thedensity=aproperty_emission[astate][distanceneg.second[argneg]].second;
        
    }
    
    
    
    return thedensity;
    
}









double logsuma(const vector<double>& logt)
{
	double expsum=0;
	int expsum_exp=0;
	const double max=maximuna(logt);
	for(size_t i=0;i<logt.size();i++)
    {
		double tmp=exp(logt[i]-max);
		for(int tt=0;tt<expsum_exp;++tt)
			tmp/=2;
        
		expsum+=tmp;
        
		while(expsum>large_num)
        {
			expsum/=2;
			++expsum_exp;
		}
	}
	return max+log(expsum)+expsum_exp*log(2);
    
}






//get the bandwidth h for the kernel estimation


double get_hvalue(const vector<double>& chosen_data)
{
    double hval;
    
    int sizen=chosen_data.size();
    
    double themean=getmean(chosen_data);
    double thevar=getvariance(chosen_data, themean);
    
    
    hval=2.0 * 1.06*pow(thevar,0.5)/pow(sizen, 0.2);
    
    
   // cout<<"my hval "<<hval;
    
    
    return hval;
    
    
    
    
}




//this function shall be used in calculating the emissions, different data size uses different h;






//here I'm going to construct the emission chart for each obs in each sequence.
//think about the data structure,

//for a single property

//these omps are written by me may be not write

vector<vector<vector<double>>>get_allstates_emission(const vector<vector<double>>&observations,const vector<vector<pair<double,double>>>& property_emission)
{
    
    const size_t len=observations.size();
    
    vector<vector<vector<double>>>allstates_emission{len};
    
    
    for (size_t i=0; i<len; i++)
    {
	allstates_emission[i].resize(observations[i].size());
        
        //cout<<i<<" critical check          "<<observations[i].size()<<endl;

         #pragma omp parallel for
        
        for (size_t p=0; p<observations[i].size(); p++)
        {
            
            vector<double>single_emi2states(2);
            
            for (size_t h=0; h<2; h++)
            {
                single_emi2states[h]=get_lnonpara(observations[i][p], h, property_emission);
                
              //  cout<<"see sigle_emi "<<single_emi2states[h]<<" "<<p<<" "<<h<<endl;
            }
            
            // #pragma omp critical
            // allstates_emission[i].push_back(single_emi4states);
            allstates_emission[i][p]=single_emi2states;
            
        }
        
        //cout<<i<<" size of each seq "<<allstates_emission[i].size()<<endl;
    }
    
    //cout<<"size of all seqs "<<allstates_emission.size()<<endl;
    
    
    return allstates_emission;
}


//since i establish this one, I need to modilfy te
/*

std::vector<std::vector<std::pair<double, double>>>get_kernel(const std::vector<std::vector<double>>&hydro_2states)
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





*/






#endif
