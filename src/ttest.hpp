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



#include <boost/math/distributions/students_t.hpp>


//in order to debug
double two_samples_t_test_equal_sd(
                                 double Sm1,
                                 double Sd1,
                                 unsigned Sn1,
                                 double Sm2,
                                 double Sd2,
                                 unsigned Sn2,
                                 double alpha)
{
    

    
    //
    // Sm1 = Sample Mean 1.
    // Sd1 = Sample Standard Deviation 1.
    // Sn1 = Sample Size 1.
    // Sm2 = Sample Mean 2.
    // Sd2 = Sample Standard Deviation 2.
    // Sn2 = Sample Size 2.
    // alpha = Significance Level.
    //
    // A Students t test applied to two sets of data.
    // We are testing the null hypothesis that the two
    // samples have the same mean and that any difference
    // if due to chance.
    // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
    //
    using namespace std;
    using namespace boost::math;
    
    
    
    double sigind;
    
    // Print header:
   // cout <<
    //"_______________________________________________\n"
    //"Student t test for two samples (equal variances)\n"
    //"_______________________________________________\n\n";
    //cout << setprecision(5);
    //cout << setw(55) << left << "Number of Observations (Sample 1)" << "=  " << Sn1 << "\n";
    //cout << setw(55) << left << "Sample 1 Mean" << "=  " << Sm1 << "\n";
    //cout << setw(55) << left << "Sample 1 Standard Deviation" << "=  " << Sd1 << "\n";
    //cout << setw(55) << left << "Number of Observations (Sample 2)" << "=  " << Sn2 << "\n";
    //cout << setw(55) << left << "Sample 2 Mean" << "=  " << Sm2 << "\n";
    //cout << setw(55) << left << "Sample 2 Standard Deviation" << "=  " << Sd2 << "\n";
    //
    // Now we can calculate and output some stats:
    //
    // Degrees of freedom:
    double v = Sn1 + Sn2 - 2;
    //cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
    // Pooled variance:
    double sp = sqrt(((Sn1-1) * Sd1 * Sd1 + (Sn2-1) * Sd2 * Sd2) / v);
    //cout << setw(55) << left << "Pooled Standard Deviation" << "=  " << sp << "\n";
    // t-statistic:
    double t_stat = (Sm1 - Sm2) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
    //cout << setw(55) << left << "T Statistic" << "=  " << t_stat << "\n";
    //
    // Define our distribution, and get the probability:
    //
    students_t dist(v);
    double q = cdf(complement(dist, fabs(t_stat)));
    //cout << setw(55) << left << "Probability that difference is due to chance" << "=  "
   // << setprecision(3) << scientific << 2 * q << "\n\n";
    //
    // Finally print out results of alternative hypothesis:
    //
//    cout << setw(55) << left <<
//    "Results for Alternative Hypothesis and alpha" << "=  "
//    << setprecision(4) << fixed << alpha << "\n\n";
//    cout << "Alternative Hypothesis              Conclusion\n";
//    cout << "Sample 1 Mean != Sample 2 Mean       " ;
//    if(q < alpha / 2){
//        cout << "NOT REJECTED\n";
//        
//    }
//    else{
//        cout << "REJECTED\n";
//        
//    }
    
    sigind=2*q;
    
//    cout << "Sample 1 Mean <  Sample 2 Mean       ";
//    if(cdf(dist, t_stat) < alpha)
//        cout << "NOT REJECTED\n";
//    else
//        cout << "REJECTED\n";
//    cout << "Sample 1 Mean >  Sample 2 Mean       ";
//    if(cdf(complement(dist, t_stat)) < alpha)
//        cout << "NOT REJECTED\n";
//    else
//        cout << "REJECTED\n";
//    cout << endl << endl;
       // cout<<sigind;
   // cout<<"I want to save for sigind"<<endl;

    return sigind;
}




//make a trick to debug




void two_samples_t_test_unequal_sd(
                                   double Sm1,
                                   double Sd1,
                                   unsigned Sn1,
                                   double Sm2,
                                   double Sd2,
                                   unsigned Sn2,
                                   double alpha)
{
    //
    // Sm1 = Sample Mean 1.
    // Sd1 = Sample Standard Deviation 1.
    // Sn1 = Sample Size 1.
    // Sm2 = Sample Mean 2.
    // Sd2 = Sample Standard Deviation 2.
    // Sn2 = Sample Size 2.
    // alpha = Significance Level.
    //
    // A Students t test applied to two sets of data.
    // We are testing the null hypothesis that the two
    // samples have the same mean and that any difference
    // if due to chance.
    // See http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm
    //
    using namespace std;
    using namespace boost::math;
    
    // Print header:
    cout <<
    "_________________________________________________\n"
    "Student t test for two samples (unequal variances)\n"
    "_________________________________________________\n\n";
    cout << setprecision(5);
    cout << setw(55) << left << "Number of Observations (Sample 1)" << "=  " << Sn1 << "\n";
    cout << setw(55) << left << "Sample 1 Mean" << "=  " << Sm1 << "\n";
    cout << setw(55) << left << "Sample 1 Standard Deviation" << "=  " << Sd1 << "\n";
    cout << setw(55) << left << "Number of Observations (Sample 2)" << "=  " << Sn2 << "\n";
    cout << setw(55) << left << "Sample 2 Mean" << "=  " << Sm2 << "\n";
    cout << setw(55) << left << "Sample 2 Standard Deviation" << "=  " << Sd2 << "\n";
    //
    // Now we can calculate and output some stats:
    //
    // Degrees of freedom:
    double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
    v *= v;
    double t1 = Sd1 * Sd1 / Sn1;
    t1 *= t1;
    t1 /=  (Sn1 - 1);
    double t2 = Sd2 * Sd2 / Sn2;
    t2 *= t2;
    t2 /= (Sn2 - 1);
    v /= (t1 + t2);
    cout << setw(55) << left << "Degrees of Freedom" << "=  " << v << "\n";
    // t-statistic:
    double t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);
    cout << setw(55) << left << "T Statistic" << "=  " << t_stat << "\n";
    //
    // Define our distribution, and get the probability:
    //
    students_t dist(v);
    double q = cdf(complement(dist, fabs(t_stat)));
    cout << setw(55) << left << "Probability that difference is due to chance" << "=  "
    << setprecision(3) << scientific << 2 * q << "\n\n";
    //
    // Finally print out results of alternative hypothesis:
    //
    cout << setw(55) << left <<
    "Results for Alternative Hypothesis and alpha" << "=  "
    << setprecision(4) << fixed << alpha << "\n\n";
    cout << "Alternative Hypothesis              Conclusion\n";
    cout << "Sample 1 Mean != Sample 2 Mean       " ;
    if(q < alpha / 2)
        cout << "NOT REJECTED\n";
    else
        cout << "REJECTED\n";
    cout << "Sample 1 Mean <  Sample 2 Mean       ";
    if(cdf(dist, t_stat) < alpha)
        cout << "NOT REJECTED\n";
    else
        cout << "REJECTED\n";
    cout << "Sample 1 Mean >  Sample 2 Mean       ";
    if(cdf(complement(dist, t_stat)) < alpha)
        cout << "NOT REJECTED\n";
    else
        cout << "REJECTED\n";
    cout << endl << endl;
}

