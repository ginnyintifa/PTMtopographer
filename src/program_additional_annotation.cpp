//
//  main.cpp
//  humansites
//
//  Created by Ginny Li on 15-7-28.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//





#include "site.hpp"
#include "aaprop.hpp"
//#include "emisson.hpp"

#include "math.hpp"
#include "input_additional_annotation.hpp"



using namespace std;


vector<int> get_domain_proportions(vector<double>pred_scores,vector<vector<string>>all_domain_bl, vector<vector<string>>all_near_domain_bl,double fdr_score)
{
    
    vector<int>domain_proportion_of_a_score(4);
    
    int all_above_score=0;
    int in_domain_above_score=0;
    int near_domain_above_score=0;
    int bt_domain_above_score=0;  //between domain means the site is not in domain and is not near domain.
    
    for (size_t i=0; i<pred_scores.size(); i++)
    {
       
        if (pred_scores[i]>=fdr_score)
        {
            all_above_score=all_above_score+1;
            
            
            if (all_domain_bl[i].size()>0)
            {
                in_domain_above_score=in_domain_above_score+1;
            }
            
            if (all_near_domain_bl[i].size()>0)
            {
                near_domain_above_score=near_domain_above_score+1;
            }
            
            if (all_domain_bl[i].size()==0&all_near_domain_bl[i].size()==0)
            {
                bt_domain_above_score=bt_domain_above_score+1;
            }
            
        }
    }
    
    
    
    domain_proportion_of_a_score[0]=in_domain_above_score;
    domain_proportion_of_a_score[1]=near_domain_above_score;
    domain_proportion_of_a_score[2]=bt_domain_above_score;
    domain_proportion_of_a_score[3]=all_above_score;
    
    
    
    cout<<in_domain_above_score<<" "<<near_domain_above_score<<" "<<bt_domain_above_score<<" "<<all_above_score<<endl;
    
    
    return domain_proportion_of_a_score;
    
}







int main(int argc, char *argv[])

{
    
    
    
    vector<string>myinput=read_input(argv[1]);
    vector<string>mypara=get_parameterlines(myinput);
    
    Para_values myparavalue=get_paravalues(mypara);
    
    
    vector<string>myproteins_id=read_tsv(myparavalue.prots_id);
    vector<vector<int>>mypositions=read_tsv_int(myparavalue.can_sites_position);
    
    
    
    
    vector<string>dom_protid=read_tsv(myparavalue.protid);
    vector<int>dom_start=read_int(myparavalue.domain_start);
    vector<int>dom_end=read_int(myparavalue.domain_end);
    vector<string>dom_name=read_tsv(myparavalue.domain_name);
    vector<string>dom_type=read_tsv(myparavalue.domain_type);
    
    
    
    vector<string> prot_cytoskeleton=read_tsv(myparavalue.prot_cytoskeleton);
    vector<string> prot_cytosol=read_tsv(myparavalue.prot_cytosol);
    vector<string> prot_endoplasmic=read_tsv(myparavalue.prot_endoplasmic);
    vector<string> prot_endosome=read_tsv(myparavalue.prot_endosome);
    vector<string> prot_extracellular=read_tsv(myparavalue.prot_extracellular);
    vector<string> prot_golgi=read_tsv(myparavalue.prot_golgi);
    vector<string> prot_lysosome=read_tsv(myparavalue.prot_lysosome);
    vector<string> prot_membrane=read_tsv(myparavalue.prot_membrane);
    vector<string> prot_mitochondrion=read_tsv(myparavalue.prot_mitochondrion);
    vector<string> prot_nucleus=read_tsv(myparavalue.prot_nucleus);
    vector<string> prot_peroxisome=read_tsv(myparavalue.prot_peroxisome);

    
    
    vector<vector<double>>rf_score=read_svm_output(myparavalue.rf_score);
    
    
    
    vector<string>head_annotation=read_tsv(myparavalue.head_annotation); //try to use read_tsv for this
    
    
    vector<string>protein_head_annotation=read_tsv(myparavalue.prot_head_annotation);
    
    //the first line of this will be the header
    
    
    
    vector<double>global_dfdr=read_tsv_double(myparavalue.global_dfdr);
    
    vector<vector<double>>prot_specific_dfdr=read_vector_double(myparavalue.prot_specific_dfdr);
    
    cout<<"got the files "<<endl;
    
    cout<<"check prot specific dfdr size "<<prot_specific_dfdr.size()<<endl;
    
    //after lunch I will program the add on dfdr part
    
    
    vector<double>all_rf_score(rf_score.size());
    
    for (size_t i=0; i<rf_score.size(); i++)
    {
        all_rf_score[i]=rf_score[i][0];
         //cout<<"predict probs "<<all_predict_probs[i]<<endl;
    }
    
    
    
    
   // vector<int>mypred_states=read_int(myparavalue.pred_site_states);

    
    cout<<myproteins_id.size()<<" look at size "<<endl;
    
    for (size_t i=0; i<myproteins_id.size(); i++)
    {
        cout<<i<<" "<<myproteins_id[i]<<endl;
    }
    
    
    vector<string>the_ids(myproteins_id.size());
    
    
    for (size_t i=0; i<myproteins_id.size(); i++)
    {
        cout<<i<<" "<<myproteins_id[i]<<endl;
        
        if (myproteins_id[i].substr(9,1)=="|")
        {
            the_ids[i]=myproteins_id[i].substr(3,6);
        }
        else
        {
            the_ids[i]=myproteins_id[i].substr(3,10);
        }
    }

    
    cout<<"getting the id "<<endl;
    
    vector<vector<vector<string>>> domain_bl{myproteins_id.size()};//for each protein and each window there may be multiple domains to which the sites belong
    
    vector<vector<vector<string>>>near_domain_bl{myproteins_id.size()};
    
    
    
    vector<vector<vector<string>>> domain_type_bl{myproteins_id.size()};
    
    vector<vector<vector<string>>> near_domain_type_bl{myproteins_id.size()};

    
    for (size_t i=0; i<myproteins_id.size(); i++)
    {
        
        vector<vector<string>>each_domain_bl{mypositions[i].size()};
        vector<vector<string>>each_domain_type_bl{mypositions[i].size()};
        vector<vector<string>>each_near_domain_bl{mypositions[i].size()};
        vector<vector<string>>each_near_domain_type_bl{mypositions[i].size()};
        
     //   cout<<"for protein "<<i<<" ******************************************************"<<endl;
        
        for (size_t d=0; d<dom_name.size(); d++)
        {
            
            size_t near_length=size_t((dom_end[d]-dom_start[d])*0.2+0.5);
            
            
            if (the_ids[i]==dom_protid[d])//check it is the same protein
            {
                for (size_t p=0; p<mypositions[i].size(); p++)
                    
                {
                    
                    
                       size_t ab_pos=mypositions[i][p]+1;//+1 because the position of domain is from 1 not 0;
                    
                        if (ab_pos>=dom_start[d]&ab_pos<=dom_end[d])
                        {
                            each_domain_bl[p].push_back(dom_name[d]);
                            each_domain_type_bl[p].push_back(dom_type[d]);
                        }
                    
                        
                        size_t star;
                        
                        if (dom_start[d]<near_length)
                        {
                            star=0;
                        }
                        
                        else
                            star=dom_start[d]-near_length;
                        
                    
                        if ((ab_pos>=star&ab_pos<dom_start[d])|(ab_pos<=(dom_end[d]+near_length)&ab_pos>dom_end[d]))
                        {
                            
                            each_near_domain_bl[p].push_back(dom_name[d]);
                            each_near_domain_type_bl[p].push_back(dom_type[d]);
                        }
                    }
                    
                }
        }
        
        
        
        domain_bl[i]=each_domain_bl;
        domain_type_bl[i]=each_domain_type_bl;
        near_domain_bl[i]=each_near_domain_bl;
        near_domain_type_bl[i]=each_near_domain_type_bl;
    }
    
    

    
    cout<<"get domain set "<<endl;
    
    
    
    
    
    
    vector<vector<string>>all_domain_bl;
    
    for (size_t i=0; i<domain_bl.size(); i++)
    {
        for (size_t p=0; p<domain_bl[i].size(); p++)
        {
                      
            all_domain_bl.push_back(domain_bl[i][p]);
        }
    }
    
    cout<<all_domain_bl.size()<<endl;
    
    vector<vector<string>>all_near_domain_bl;
    
    for (size_t i=0; i<near_domain_bl.size(); i++)
    {
        for (size_t p=0; p<near_domain_bl[i].size(); p++)
        {
            
            all_near_domain_bl.push_back(near_domain_bl[i][p]);
        }
    }
    
    cout<<all_near_domain_bl.size()<<endl;
    
   
    
    
    
    ofstream outfile;
    
    vector<vector<string>>each_protein_localization{myproteins_id.size()};
    
    
    for (size_t i=0; i<myproteins_id.size(); i++)
    {
        for (size_t h=0; h<prot_cytoskeleton.size(); h++)
        {
            if (the_ids[i]==prot_cytoskeleton[h])
            {
                each_protein_localization[i].push_back("cytoskeleton");
            }
        }
        
        for (size_t h=0; h<prot_cytosol.size(); h++)
        {
            if (the_ids[i]==prot_cytosol[h])
            {
                each_protein_localization[i].push_back("cytosol");
            }
        }
        
        for (size_t h=0; h<prot_endoplasmic.size(); h++)
        {
            if (the_ids[i]==prot_endoplasmic[h])
            {
                each_protein_localization[i].push_back("endoplasmic");
            }
        }
        for (size_t h=0; h<prot_endosome.size(); h++)
        {
            if (the_ids[i]==prot_endosome[h])
            {
                each_protein_localization[i].push_back("endosome");
            }
        }
        
        for (size_t h=0; h<prot_extracellular.size(); h++)
        {
            if (the_ids[i]==prot_extracellular[h])
            {
                each_protein_localization[i].push_back("extracellular");
            }
        }
        for (size_t h=0; h<prot_golgi.size(); h++)
        {
            if (the_ids[i]==prot_golgi[h])
            {
                each_protein_localization[i].push_back("golgi");
            }
        }
        
        for (size_t h=0; h<prot_lysosome.size(); h++)
        {
            if (the_ids[i]==prot_lysosome[h])
            {
                each_protein_localization[i].push_back("lysosome");
            }
        }
        
        for (size_t h=0; h<prot_membrane.size(); h++)
        {
            if (the_ids[i]==prot_membrane[h])
            {
                each_protein_localization[i].push_back("membrane");
            }
        }
        
        for (size_t h=0; h<prot_mitochondrion.size(); h++)
        {
            if (the_ids[i]==prot_mitochondrion[h])
            {
                each_protein_localization[i].push_back("mitochondrion");
            }
        }
        for (size_t h=0; h<prot_nucleus.size(); h++)
        {
            if (the_ids[i]==prot_nucleus[h])
            {
                each_protein_localization[i].push_back("nucleus");
            }
        }
        
        for (size_t h=0; h<prot_peroxisome.size(); h++)
        {
            if (the_ids[i]==prot_peroxisome[h])
            {
                each_protein_localization[i].push_back("peroxisome");
            }
        }
        
    }
    
    
    
    vector<int>site_prot;
    
    for (size_t i=0; i<mypositions.size(); i++)
    {
        for (size_t h=0; h<mypositions[i].size(); h++)
        {
            site_prot.push_back(i);
        }
    }
    
    cout<<"NEW VECTOR SIZE "<<all_rf_score.size()<<" "<<site_prot.size()<<endl;
    
    
    
    
    
    outfile.open("site_annotation.tsv");
    
    
    outfile<<head_annotation[0];
    outfile<<"prediction_score"<<"\t"<<"global_dfdr"<<"\t"<<"pspc_dfdr"<<"\t"<<"domain_in"<<"\t"<<"domain_near"<<endl;
    
    for (size_t i=0; i<head_annotation.size()-1; i++)
    {
       
           outfile<<head_annotation[i+1]<<all_rf_score[i]<<"\t";
        
        int which_dfdr=all_rf_score[i]*1000-1;
        
        
        if (which_dfdr<0)
        {
            outfile<<1<<"\t"<<1<<"\t";
        }
        else
        {
            outfile<<(global_dfdr[which_dfdr]>1?1:global_dfdr[which_dfdr])<<"\t"<<(prot_specific_dfdr[site_prot[i]][which_dfdr]>1?1:prot_specific_dfdr[site_prot[i]][which_dfdr])<<"\t";
            
        }
        
        
            for (size_t h=0; h<all_domain_bl[i].size(); h++)
            {
                outfile<<all_domain_bl[i][h]<<" ";
            }
            outfile<<"\t";
            
            for (size_t h=0; h<all_near_domain_bl[i].size(); h++)
            {
                outfile<<all_near_domain_bl[i][h]<<" ";
            }
        
            outfile<<endl;
        
    }
    
    outfile.close();
    
    
    
    
    outfile.open("protein_annotation.tsv");
    
    
    outfile<<protein_head_annotation[0];
    outfile<<"subcellular_location"<<endl;
    
    cout<<head_annotation.size()<<endl;
    
    for (size_t i=0; i<protein_head_annotation.size()-1; i++)
    {
       // cout<<i<<" not wrong "<<endl;
        outfile<<protein_head_annotation[i+1];
        
        for (size_t h=0; h<each_protein_localization[i].size(); h++)
        {
            outfile<<each_protein_localization[i][h]<<" ";
        }
        
        outfile<<endl;
       
    }
    
    outfile.close();
    
    
    vector<int> proportion_first=get_domain_proportions(all_rf_score, all_domain_bl,all_near_domain_bl,myparavalue.first_fdr_score);
    
    
    
    vector<int> proportion_second=get_domain_proportions(all_rf_score, all_domain_bl,all_near_domain_bl,myparavalue.second_fdr_score);
    
    
    vector<int> proportion_third=get_domain_proportions(all_rf_score, all_domain_bl,all_near_domain_bl,myparavalue.third_fdr_score);
    
    
    
    
    vector<int> proportion_all=get_domain_proportions(all_rf_score, all_domain_bl,all_near_domain_bl,0);
    
    
    
}

