//
//  aaprop.hpp
//  humansites
//
//  Created by Ginny Li on 15-7-29.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//

#ifndef aaprop_hpp
#define aaprop_hpp

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


//now I have to add in pka1 and pka2

//change this to aa_properties

struct aa_properties{
    char aaname;
    
    double isoelec;
    double pka1;
    double pka2;
    double hydro;
    double helixpro;
    double stericpar;
    double polar;
    double volume;
    double sheetpro;
    double aafreq;
    
};


struct mean_properties{

    
    double mean_pka1;
    double mean_hydro;
    double mean_helixpro;

    double mean_stericpar;
    double mean_polar;
    double mean_volume;
    double mean_sheetpro;
 

    
};




/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////get the 3 properties for each block, free of sites information

//A codebook for each amino acid

//Xchange

//modify this part in the afternoon

//add the properties in

std::vector<aa_properties> aminoacid_properties(24);

void assign_properties(std::vector<aa_properties>& aminoacid_properties)
{
    aminoacid_properties[23].aaname='Z';
    aminoacid_properties[23].isoelec=0;
    aminoacid_properties[23].hydro=0;
    aminoacid_properties[23].pka1=0;
    aminoacid_properties[23].pka2=0;
    aminoacid_properties[23].helixpro=0;
    aminoacid_properties[23].stericpar=0;
    aminoacid_properties[23].polar=0;
    aminoacid_properties[23].volume=0;
    aminoacid_properties[23].sheetpro=0;
    
    aminoacid_properties[22].aaname='U';
    aminoacid_properties[22].isoelec=0;
    aminoacid_properties[22].hydro=0;
    aminoacid_properties[22].pka1=0;
    aminoacid_properties[22].pka2=0;
    aminoacid_properties[22].helixpro=0;
    aminoacid_properties[22].stericpar=0;
    aminoacid_properties[22].polar=0;
    aminoacid_properties[22].volume=0;
    aminoacid_properties[22].sheetpro=0;
    
    
    
    aminoacid_properties[21].aaname='B';
    aminoacid_properties[21].isoelec=0;
    aminoacid_properties[21].hydro=0;
    aminoacid_properties[21].pka1=0;
    aminoacid_properties[21].pka2=0;
    aminoacid_properties[21].helixpro=0;
    aminoacid_properties[21].stericpar=0;
    aminoacid_properties[21].polar=0;
    aminoacid_properties[21].volume=0;
    aminoacid_properties[21].sheetpro=0;
    
    
    
    
    aminoacid_properties[20].aaname='X';
    aminoacid_properties[20].isoelec=0;
    aminoacid_properties[20].hydro=0;
    aminoacid_properties[20].pka1=0;
    aminoacid_properties[20].pka2=0;
    aminoacid_properties[20].helixpro=0;
    aminoacid_properties[20].stericpar=0;
    aminoacid_properties[20].polar=0;
    aminoacid_properties[20].volume=0;
    aminoacid_properties[20].sheetpro=0;
    
    
    
    aminoacid_properties[0].aaname='G';
    aminoacid_properties[0].isoelec=5.97;
    aminoacid_properties[0].hydro=0.00;
    aminoacid_properties[0].pka1=2.34;
    aminoacid_properties[0].pka2=9.60;
    aminoacid_properties[0].helixpro=0.56;
    aminoacid_properties[0].stericpar=0;
    aminoacid_properties[0].polar=0;
    aminoacid_properties[0].volume=0;
    aminoacid_properties[0].sheetpro=0.92;
    
    
    
    aminoacid_properties[1].aaname='A';
    aminoacid_properties[1].isoelec=6.01;
    aminoacid_properties[1].hydro=0.31;
    aminoacid_properties[1].pka1=2.34;
    aminoacid_properties[1].pka2=9.69;
    aminoacid_properties[1].helixpro=1.29;
    aminoacid_properties[1].stericpar=1.28;
    aminoacid_properties[1].polar=0.05;
    aminoacid_properties[1].volume=1.00;
    aminoacid_properties[1].sheetpro=0.90;
    
    
    aminoacid_properties[2].aaname='L';
    aminoacid_properties[2].isoelec=5.98;
    aminoacid_properties[2].hydro=1.70;
    aminoacid_properties[2].pka1=2.36;
    aminoacid_properties[2].pka2=9.60;
    aminoacid_properties[2].helixpro=1.30;
    aminoacid_properties[2].stericpar=2.59;
    aminoacid_properties[2].polar=0.19;
    aminoacid_properties[2].volume=4.00;
    aminoacid_properties[2].sheetpro=1.02;
    
    
    
    aminoacid_properties[3].aaname='V';
    aminoacid_properties[3].isoelec=5.97;
    aminoacid_properties[3].hydro=1.22;
    aminoacid_properties[3].pka1=2.32;
    aminoacid_properties[3].pka2=9.62;
    aminoacid_properties[3].helixpro=0.91;
    aminoacid_properties[3].stericpar=3.67;
    aminoacid_properties[3].polar=0.14;
    aminoacid_properties[3].volume=3.00;
    aminoacid_properties[3].sheetpro=1.49;
    
    
    
    aminoacid_properties[4].aaname='I';
    aminoacid_properties[4].isoelec=6.02;
    aminoacid_properties[4].hydro=1.80;
    aminoacid_properties[4].pka1=2.36;
    aminoacid_properties[4].pka2=9.60;
    aminoacid_properties[4].helixpro=0.97;
    aminoacid_properties[4].stericpar=4.19;
    aminoacid_properties[4].polar=0.19;
    aminoacid_properties[4].volume=4.00;
    aminoacid_properties[4].sheetpro=1.45;
    
    
    
    aminoacid_properties[5].aaname='F';
    aminoacid_properties[5].isoelec=5.48;
    aminoacid_properties[5].hydro=1.79;
    aminoacid_properties[5].pka1=1.83;
    aminoacid_properties[5].pka2=9.13;
    aminoacid_properties[5].helixpro=0.30;
    aminoacid_properties[5].stericpar=2.94;
    aminoacid_properties[5].polar=0.29;
    aminoacid_properties[5].volume=5.89;
    aminoacid_properties[5].sheetpro=0.38;
    
    
    
    aminoacid_properties[6].aaname='W';
    aminoacid_properties[6].isoelec=5.89;
    aminoacid_properties[6].hydro=2.25;
    aminoacid_properties[6].pka1=2.83;
    aminoacid_properties[6].pka2=9.39;
    aminoacid_properties[6].helixpro=0.99;
    aminoacid_properties[6].stericpar=3.21;
    aminoacid_properties[6].polar=0.41;
    aminoacid_properties[6].volume=8.08;
    aminoacid_properties[6].sheetpro=1.14;
    
    
    
    aminoacid_properties[7].aaname='Y';
    aminoacid_properties[7].isoelec=5.67;
    aminoacid_properties[7].hydro=0.96;
    aminoacid_properties[7].pka1=2.20;
    aminoacid_properties[7].pka2=9.11;
    aminoacid_properties[7].helixpro=0.72;
    aminoacid_properties[7].stericpar=2.94;
    aminoacid_properties[7].polar=0.30;
    aminoacid_properties[7].volume=6.47;
    aminoacid_properties[7].sheetpro=1.25;
    
    
    
    aminoacid_properties[8].aaname='P';
    aminoacid_properties[8].isoelec=6.48;
    aminoacid_properties[8].hydro=0.72;
    aminoacid_properties[8].pka1=1.99;
    aminoacid_properties[8].pka2=10.60;
    aminoacid_properties[8].helixpro=0.52;
    aminoacid_properties[8].stericpar=2.67;
    aminoacid_properties[8].polar=0.00;
    aminoacid_properties[8].volume=2.72;
    aminoacid_properties[8].sheetpro=0.64;
    
    
    
    aminoacid_properties[9].aaname='M';
    aminoacid_properties[9].isoelec=5.47;
    aminoacid_properties[9].hydro=1.23;
    aminoacid_properties[9].pka1=2.28;
    aminoacid_properties[9].pka2=9.21;
    aminoacid_properties[9].helixpro=1.47;
    aminoacid_properties[9].stericpar=2.35;
    aminoacid_properties[9].polar=0.22;
    aminoacid_properties[9].volume=4.43;
    aminoacid_properties[9].sheetpro=0.97;
    
    
    
    aminoacid_properties[10].aaname='C';
    aminoacid_properties[10].isoelec=5.07;
    aminoacid_properties[10].hydro=1.54;
    aminoacid_properties[10].pka1=1.96;
    aminoacid_properties[10].pka2=8.18;
    aminoacid_properties[10].helixpro=1.11;
    aminoacid_properties[10].stericpar=1.77;
    aminoacid_properties[10].polar=0.13;
    aminoacid_properties[10].volume=2.43;
    aminoacid_properties[10].sheetpro=0.74;
    
    
    
    aminoacid_properties[11].aaname='T';
    aminoacid_properties[11].isoelec=5.87;
    aminoacid_properties[11].hydro=0.26;
    aminoacid_properties[11].pka1=2.09;
    aminoacid_properties[11].pka2=9.10;
    aminoacid_properties[11].helixpro=0.82;
    aminoacid_properties[11].stericpar=3.03;
    aminoacid_properties[11].polar=0.11;
    aminoacid_properties[11].volume=2.60;
    aminoacid_properties[11].sheetpro=1.21;
    
    
    
    aminoacid_properties[12].aaname='S';
    aminoacid_properties[12].isoelec=5.68;
    aminoacid_properties[12].hydro=-0.04;
    aminoacid_properties[12].pka1=2.21;
    aminoacid_properties[12].pka2=9.15;
    aminoacid_properties[12].helixpro=0.82;
    aminoacid_properties[12].stericpar=1.31;
    aminoacid_properties[12].polar=0.06;
    aminoacid_properties[12].volume=1.60;
    aminoacid_properties[12].sheetpro=0.95;
    
    
    
    aminoacid_properties[13].aaname='N';
    aminoacid_properties[13].isoelec=5.41;
    aminoacid_properties[13].hydro=-0.60;
    aminoacid_properties[13].pka1=2.02;
    aminoacid_properties[13].pka2=8.80;
    aminoacid_properties[13].helixpro=0.90;
    aminoacid_properties[13].stericpar=1.60;
    aminoacid_properties[13].polar=0.13;
    aminoacid_properties[13].volume=2.95;
    aminoacid_properties[13].sheetpro=0.76;
    
    
    
    aminoacid_properties[14].aaname='Q';
    aminoacid_properties[14].isoelec=5.65;
    aminoacid_properties[14].hydro=-0.22;
    aminoacid_properties[14].pka1=2.17;
    aminoacid_properties[14].pka2=9.13;
    aminoacid_properties[14].helixpro=1.27;
    aminoacid_properties[14].stericpar=1.56;
    aminoacid_properties[14].polar=0.18;
    aminoacid_properties[14].volume=3.95;
    aminoacid_properties[14].sheetpro=0.80;
    
    
    
    aminoacid_properties[15].aaname='D';
    aminoacid_properties[15].isoelec=2.77;
    aminoacid_properties[15].hydro=-0.77;
    aminoacid_properties[15].pka1=1.88;
    aminoacid_properties[15].pka2=9.60;
    aminoacid_properties[15].helixpro=1.04;
    aminoacid_properties[15].stericpar=1.60;
    aminoacid_properties[15].polar=0.11;
    aminoacid_properties[15].volume=2.78;
    aminoacid_properties[15].sheetpro=0.72;
    
    
    
    aminoacid_properties[16].aaname='E';
    aminoacid_properties[16].isoelec=3.22;
    aminoacid_properties[16].hydro=-0.64;
    aminoacid_properties[16].pka1=2.19;
    aminoacid_properties[16].pka2=9.67;
    aminoacid_properties[16].helixpro=1.44;
    aminoacid_properties[16].stericpar=1.56;
    aminoacid_properties[16].polar=0.15;
    aminoacid_properties[16].volume=3.78;
    aminoacid_properties[16].sheetpro=0.75;
    
    
    
    aminoacid_properties[17].aaname='H';
    aminoacid_properties[17].isoelec=7.59;
    aminoacid_properties[17].hydro=0.13;
    aminoacid_properties[17].pka1=1.82;
    aminoacid_properties[17].pka2=9.17;
    aminoacid_properties[17].helixpro=1.22;
    aminoacid_properties[17].stericpar=2.99;
    aminoacid_properties[17].polar=0.23;
    aminoacid_properties[17].volume=4.66;
    aminoacid_properties[17].sheetpro=1.08;
    
    
    
    aminoacid_properties[18].aaname='K';
    aminoacid_properties[18].isoelec=9.74;
    aminoacid_properties[18].hydro=-0.99;
    aminoacid_properties[18].pka1=2.18;
    aminoacid_properties[18].pka2=8.95;
    aminoacid_properties[18].helixpro=1.23;
    aminoacid_properties[18].stericpar=1.89;
    aminoacid_properties[18].polar=0.22;
    aminoacid_properties[18].volume=4.77;
    aminoacid_properties[18].sheetpro=0.77;
    
    
    
    aminoacid_properties[19].aaname='R';
    aminoacid_properties[19].isoelec=10.76;
    aminoacid_properties[19].hydro=-1.01;
    aminoacid_properties[19].pka1=2.17;
    aminoacid_properties[19].pka2=9.04;
    aminoacid_properties[19].helixpro=0.96;
    aminoacid_properties[19].stericpar=2.34;
    aminoacid_properties[19].polar=0.29;
    aminoacid_properties[19].volume=6.13;
    aminoacid_properties[19].sheetpro=0.99;
    
    
}




//A function to check which site corresponds to which number

std::vector<int> get_aawhich(std::string pred_site)
{
    using namespace std;
    
    vector<int>all_aawhich;
    
    for (size_t i=0; i<pred_site.size(); i++)
    {
        
    int aawhich;
    
    switch (pred_site[i])
    {
        case 'z':
            aawhich=23;
        case 'u':
            aawhich=22;
        case 'b':
            aawhich=21;
        case 'x':
            aawhich=20;
        case 'g':
            aawhich=0;
            break;
        case 'a':
            aawhich=1;
            break;
        case 'l':
            aawhich=2;
            break;
        case 'v':
            aawhich=3;
            break;
        case 'i':
            aawhich=4;
            break;
            
        case 'f':
            aawhich=5;
            break;
            
        case 'w':
            aawhich=6;
            break;
        case 'y':
            aawhich=7;
            break;
        case 'p':
            aawhich=8;
            break;
            
        case 'm':
            aawhich=9;
            break;
        case 'c':
            aawhich=10;
            break;
            
        case 't':
            aawhich=11;
            break;
            
        case 's':
            aawhich=12;
            break;
            
        case 'n':
            aawhich=13;
            break;
            
        case 'q':
            aawhich=14;
            break;
            
        case 'd':
            aawhich=15;
            break;
            
        case 'e':
            aawhich=16;
            break;
        case 'h':
            aawhich=17;
            break;
        case 'k':
            aawhich=18;
            break;
            
        case 'r':
            aawhich=19;
            break;
            
        default:
            cout<<"invalid aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa "<<pred_site[i]<<endl;
        
    }
        all_aawhich.push_back(aawhich);
        
    }
    
    return all_aawhich;
}






std::vector<int>grab_aa_vector(const std::vector<aa_properties> &aminoacid_properties,const std::string& awindow)
{
    using namespace std;
    vector<int> aa_vector(24);
    for (size_t k=0; k<24; k++) {
        aa_vector[k]=0;
    }
    
    
    for (int i=0; i<awindow.size(); i++)
    {
        
        char myaa=awindow[i];
        
        switch (myaa)
        {
            case 'Z':
                aa_vector[23]+=1;
                
            case 'U':
                aa_vector[22]+=1;
                
            case 'B':
                aa_vector[21]+=1;
            case 'X':
                aa_vector[20]+=1;
            case 'G':
                aa_vector[0]+=1;
                break;
            case 'A':
                aa_vector[1]+=1;
                break;
            case 'L':
                aa_vector[2]+=1;
                break;
            case 'V':
                aa_vector[3]+=1;
                break;
                
            case 'I':
                aa_vector[4]+=1;
                break;
                
            case 'F':
                aa_vector[5]+=1;
                break;
                
            case 'W':
                aa_vector[6]+=1;
                break;
            case 'Y':
                aa_vector[7]+=1;
                break;
            case 'P':
                aa_vector[8]+=1;
                break;
                
            case 'M':
                aa_vector[9]+=1;
                break;
            case 'C':
                aa_vector[10]+=1;
                break;
                
            case 'T':
                aa_vector[11]+=1;
                break;
                
            case 'S':
                aa_vector[12]+=1;
                break;
                
            case 'N':
                aa_vector[13]+=1;
                break;
                
            case 'Q':
                aa_vector[14]+=1;
                break;
                
            case 'D':
                aa_vector[15]+=1;
                break;
                
            case 'E':
                aa_vector[16]+=1;
                break;
            case 'H':
                aa_vector[17]+=1;
                break;
            case 'K':
                aa_vector[18]+=1;
                break;
                
            case 'R':
                aa_vector[19]+=1;
                break;
                
            default:
                cout<<"invalid aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa "<<myaa<<endl;
                
        }
        
        
    }
    
    return aa_vector;
}

//make the following function a single one---in train emisson program

//Xchange
//double get_average_hyd(const std::vector<aa_iso_hyd_helix> &aminoacid_properties, const std::vector<int>& aa_vector)
//{
//    using namespace std;
//    
//    
//    double average_hyd;
//    double sum_num=0;
//    
//    
//    
//    double sum_hydrophobicity=0;
//    
//    
//    for (int i=0;i<20;i++)
//    {
//        sum_hydrophobicity+=aa_vector[i]*aminoacid_properties[i].hydro;
//        
//        sum_num+=aa_vector[i];
//        
//    }
//     average_hyd=sum_hydrophobicity/sum_num;
//    
//    
//    return average_hyd;
//    
//}
//
//Xchange


std::vector<double>each_hydro;
std::vector<double>each_pka1;
std::vector<double>each_helixpro;
std::vector<double>each_stericpar;
std::vector<double>each_polar;
std::vector<double>each_volume;
std::vector<double>each_sheetpro;

//not sure where should I assign values to them
void assign_each_property(std::vector<aa_properties>aminoacid_properties,
                          std::vector<double>& each_hydro,
                          std::vector<double>& each_pka1,
                          std::vector<double>& each_helixpro,
                          std::vector<double>& each_stericpar,
                          std::vector<double>& each_polar,
                          std::vector<double>& each_volume,
                          std::vector<double>& each_sheetpro)

{
    for(size_t i=0;i<aminoacid_properties.size();i++)
    {
        each_hydro.push_back(aminoacid_properties[i].hydro);
        each_pka1.push_back(aminoacid_properties[i].pka1);
        each_helixpro.push_back(aminoacid_properties[i].helixpro);
        each_stericpar.push_back(aminoacid_properties[i].stericpar);
        each_polar.push_back(aminoacid_properties[i].polar);
        each_volume.push_back(aminoacid_properties[i].volume);
        each_sheetpro.push_back(aminoacid_properties[i].sheetpro);
    }

}


//notice the significant change in the input parameter


 double get_average_property(const std::vector<double> &each_property, const std::vector<int>& aa_vector)
{
    using namespace std;
    
    
    double average_property;
    
    double sum_property=0;
    double sum_num=0;
    
    for (int i=0;i<20;i++)
    {
        sum_property+=aa_vector[i]*each_property[i];
        
        
        sum_num+=aa_vector[i];
    }
    
    
     average_property=sum_property/sum_num;
    
    return average_property;
    
}




//Xchange

std::vector<std::vector<std::vector<int>>>get_allseqs__aavector(const std::vector<std::vector<std::string>>&window_for_all, const std::vector<aa_properties> & aminoacid_properties)
{
    using namespace std;
    
    vector<vector<vector<int>>> aa_compo_forall;
    
    
    for (size_t i=0; i<window_for_all.size(); i++)
    {
        vector<vector<int>>aa_vector_aseq;
        for (size_t p=0; p<window_for_all[i].size(); p++)
        {
        
            vector<int>aa_vector=grab_aa_vector(aminoacid_properties,window_for_all[i][p]);
            aa_vector_aseq.push_back(aa_vector);
        }
        
        aa_compo_forall.push_back(aa_vector_aseq);
    }
    
    
    return aa_compo_forall;
    

    
    
}



//get hydro and isoelectric for all seqs

//Xchange

std::vector<std::vector<double>>get_allseqs_property(const std::vector<std::vector<std::vector<int>>>& aacompoforall,const std::vector<double>  & each_property ,const std::vector<double>& property_means)
{
    using namespace std;
    
    vector<vector<double>> property_forall;
    
    
    for (size_t i=0; i<aacompoforall.size(); i++)
    {
        vector<double>property_aseq;
        for (size_t p=0; p<aacompoforall[i].size(); p++)
        {
            double property_awindow=get_average_property(each_property, aacompoforall[i][p])-property_means[i];
                   property_aseq.push_back(property_awindow);
        }
                property_forall.push_back(property_aseq);
    }
    
    
    return property_forall;

    
}
//
//std::vector<std::vector<double>>get_allseqs_pka1(const std::vector<std::vector<std::vector<int>>>& aacompoforall,const std::vector<aa_iso_hyd_helix> & aminoacid_properties,const std::vector<double>& pka1_means)
//{
//    using namespace std;
//    
//    vector<vector<double>> pka1_forall;
//    
//    
//    for (size_t i=0; i<aacompoforall.size(); i++)
//    {
//        vector<double>pka1_aseq;
//        for (size_t p=0; p<aacompoforall[i].size(); p++)
//        {
//            double pka1_awindow=get_average_pka1(aminoacid_properties, aacompoforall[i][p])-pka1_means[i];
//                pka1_aseq.push_back(pka1_awindow);
//        }
//        
//        pka1_forall.push_back(pka1_aseq);
//    }
//    
//    
//    return pka1_forall;
//    
//    
//}
//
//
//
//
//std::vector<std::vector<double>>get_allseqs_helix(const std::vector<std::vector<std::vector<int>>>& aacompoforall,const std::vector<aa_iso_hyd_helix> & aminoacid_properties,const std::vector<double>& helix_means)
//{
//    using namespace std;
//    
//    vector<vector<double>> helix_forall;
//    
//    
//    for (size_t i=0; i<aacompoforall.size(); i++)
//    {
//        vector<double>helix_aseq;
//        for (size_t p=0; p<aacompoforall[i].size(); p++)
//        {
//            double helix_awindow=get_average_helix(aminoacid_properties, aacompoforall[i][p])-helix_means[i];
//            helix_aseq.push_back(helix_awindow);
//        }
//        
//        helix_forall.push_back(helix_aseq);
//    }
//    
//    
//    return helix_forall;
//    
//    
//}
//
//
//











//ok what I do now is to use the aa compo for all and siteswith windows to program so that I will get for each protein, for each amino acid

//so the data structure is

//vector<vector<int>>-----protein-aa


std::vector<std::vector<std::pair<int,int>>> get_windows_ratio(std::vector<std::vector<std::vector<int>>> aacompoforall, std::vector<std::pair<prot_type,std::vector<int>>>allseqs_blockwithsites)
{
    using namespace std;
    
    vector<vector<pair<int,int>>> windowsratio{aacompoforall.size()};
    
    
    
    
    for (size_t i=0;i<aacompoforall.size();i++)
    {
        vector<pair<int,int>>aacounts(20);
        
        for (size_t q=0; q<20; q++)
        {
            aacounts[q].first=0;
        }
        
        
        
        for (size_t p=0; p<allseqs_blockwithsites[i].second.size(); p++)
        {
            for (size_t h=0; h<20; h++)
            {
                if(aacompoforall[i][allseqs_blockwithsites[i].second[p]][h]>0)
                aacounts[h].first=aacounts[h].first+1;
            }
            
        }
        
        for (size_t j=0;j<20;j++)
        {
            aacounts[j].second=allseqs_blockwithsites[i].second.size()-aacounts[j].first;
        }
        
        windowsratio[i]=aacounts;
    }
    
    
    return windowsratio;
    
    
    
}

std::vector<std::vector<std::pair<int, int>>>get_statesrates(std::vector<std::vector<std::vector<int>>> aacompoforall,std::vector<std::vector<int>> windowstates)
{
    using namespace std;
    
    vector<vector<pair<int,int>>>states_windowrates{2};
    vector<int> statesall(2);
    
        for (size_t k=0;k<2;k++)
        {
            
            vector<pair<int,int>>aacounts(20);
            
            for (size_t q=0; q<20; q++)
            {
                aacounts[q].first=0;
                aacounts[q].second=0;
            }
            
            statesall[k]=0;
            
            for (size_t i=0;i<aacompoforall.size();i++)
            {
              for (size_t p=0; p<aacompoforall[i].size(); p++)
               {
                 if (windowstates[i][p]==k)
                {
                  for (size_t h=0; h<20; h++)
                  {
                   aacounts[h].first=aacounts[h].first+aacompoforall[i][p][h];
                   
                   statesall[k]=statesall[k]+aacompoforall[i][p][h];
   
                  }
                }
            
               }

         }
            
            for (size_t e=0; e<20; e++)
            {
                aacounts[e].second=statesall[k];
            }
            
            states_windowrates[k]=aacounts;
       }
    

    return states_windowrates;
    
    
}


//notice the change of input
std::vector<mean_properties> get_means_all(const std::vector<double> & each_hydro, const std::vector<double> & each_pka1,const std::vector<double> & each_helixpro,const std::vector<double> & each_stericpar,const std::vector<double> & each_polar,const std::vector<double> & each_volume,const std::vector<double> & each_sheetpro, const std::vector<std::vector<int>>& aa_vectorforall)
{
    using namespace std;
    
    vector<mean_properties> seq_means(aa_vectorforall.size());
    
    
    for (size_t i=0;i<aa_vectorforall.size();i++)
    {
        seq_means[i].mean_hydro=get_average_property(each_hydro, aa_vectorforall[i]);
        
        seq_means[i].mean_pka1=get_average_property(each_pka1, aa_vectorforall[i]);
        
        seq_means[i].mean_helixpro=get_average_property(each_helixpro, aa_vectorforall[i]);
        
        seq_means[i].mean_stericpar=get_average_property(each_stericpar, aa_vectorforall[i]);
        
        seq_means[i].mean_polar=get_average_property(each_polar, aa_vectorforall[i]);
        
        seq_means[i].mean_volume=get_average_property(each_volume, aa_vectorforall[i]);
        
        seq_means[i].mean_sheetpro=get_average_property(each_sheetpro, aa_vectorforall[i]);
        


        
        
    }
    
    
    
    
    return seq_means;
    
    
    
    
}




std::vector<double> get_means(const std::vector<double> & each_property,const std::vector<std::vector<int>>& aa_vectorforall)
{
    using namespace std;
    
    vector<double> seq_means(aa_vectorforall.size());
    
    
    for (size_t i=0;i<aa_vectorforall.size();i++)
    {
        seq_means[i]=get_average_property(each_property, aa_vectorforall[i]);
        
        
    }
    
    return seq_means;
    
    
}



//this one maybe used multiple times in main.cpp


std::vector<std::vector<double>>get_properties_2states(const std::vector<std::vector<double>>&hyd_forall,const std::vector<std::vector<int>>&windows_states)//input make it a simple vector)

{
    using namespace std;
    vector<vector<double>>myproperty_2states{2};
    
    for (size_t k=0; k<2; k++)
    {
      for (size_t i=0; i<hyd_forall.size(); i++)
      {
    
    
      for (size_t p=0; p<hyd_forall[i].size(); p++)
      {
        
        if (windows_states[i][p]==k)
        {
        myproperty_2states[k].push_back(hyd_forall[i][p]);//data structure.
        }
            
      }
            
            
      }
    }
    

    return myproperty_2states;
}




#endif
