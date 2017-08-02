//
//  main.cpp
//  humansites
//
//  Created by Ginny Li on 15-7-28.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//



#include "site.hpp"
#include "aaprop.hpp"
#include "emisson.hpp"
#include "math.hpp"
#include "input_prediction_summary.hpp"
#include "ttest.hpp"
#include "chitest.hpp"


using namespace std;



vector<int>get_tf_table(vector<double>global_dfdr, double given_dfdr,vector<double>predicted_score, vector<int>known_state)
{
    vector<double> abs_dif(global_dfdr.size());
    
    
    for (size_t i=0; i<global_dfdr.size(); i++)
    {
        abs_dif[i]=abs(given_dfdr-global_dfdr[i]);
    }
    
    int which_arg=get_argminimuma(abs_dif);
    
    
    double which_score=double(which_arg)/double(10000);
    
    //cout<<which_score<<endl;
    
    
    vector<int>tf_table(4);
    
    tf_table[0]=tf_table[1]=tf_table[2]=tf_table[3]=0;
    
    for (size_t i=0; i<predicted_score.size(); i++)
    {
        if (predicted_score[i]>=which_score)
        {
            if (known_state[i]==1)
            {
                tf_table[0]=tf_table[0]+1;
            }
            else
            {
                tf_table[1]=tf_table[1]+1;
            }
        }
        else
        {
            if (known_state[i]==1)
            {
                tf_table[2]=tf_table[2]+1;
            }
            else
            {
                tf_table[3]=tf_table[3]+1;
            }
        }
    }
    
    
    
    return tf_table;
    
    
    
}


vector<int>get_tf_table_via_score(double score, vector<double>predicted_score, vector<int>known_state)
{
    
    
    double which_score=score;
   // cout<<which_score<<endl;
    
    
    vector<int>tf_table(4);
    
    tf_table[0]=tf_table[1]=tf_table[2]=tf_table[3]=0;
    
    for (size_t i=0; i<predicted_score.size(); i++)
    {
        if (predicted_score[i]>=which_score)
        {
            if (known_state[i]==1)
            {
                tf_table[0]=tf_table[0]+1;
            }
            else
            {
                tf_table[1]=tf_table[1]+1;
            }
        }
        else
        {
            if (known_state[i]==1)
            {
                tf_table[2]=tf_table[2]+1;
            }
            else
            {
                tf_table[3]=tf_table[3]+1;
            }
        }
    }
    
    
    
    return tf_table;
    
    
    
}






double get_score_for_max_or(vector<double>predicted_score,vector<int>known_state)
{
    vector<double>or_cal(999);
    
    
    for (size_t h=1; h<1000; h++)
    {
        
        double which_score=double(h)/double(1000);
        
        vector<int>tf_table(4);
        
        tf_table[0]=tf_table[1]=tf_table[2]=tf_table[3]=0;
        for (size_t i=0; i<predicted_score.size(); i++)
        {
            if (predicted_score[i]>=which_score)
            {
                if (known_state[i]==1)
                {
                    tf_table[0]=tf_table[0]+1;
                }
                else
                {
                    tf_table[1]=tf_table[1]+1;
                }
            }
            else
            {
                if (known_state[i]==1)
                {
                    tf_table[2]=tf_table[2]+1;
                }
                else
                {
                    tf_table[3]=tf_table[3]+1;
                }
            }
            
        
        }
        
        or_cal[h-1]=double(tf_table[0])/double(tf_table[1])*(double(tf_table[3])/double(tf_table[2]));


    }
    
        int arg_score=get_argmaximuma(or_cal);
    
    
   // cout<<arg_score<<" "<<or_cal[arg_score]<<endl;
    
    return double(arg_score+1)/double(10000);
    
    
    
}







int main(int argc, char *argv[])

{
    
    cout<<"Start to run"<<endl;
    
    
    
    vector<string>myinput=read_input(argv[1]);
    
    vector<string>mypara=get_parameterlines(myinput);
    
    Para_values myparavalue=get_paravalues(mypara);
    
    cout<<"Get my paravalue "<<endl;
    
    
    
    
    
    
    cout<<myparavalue.pred_probs<<endl;
    
    
    
    vector<vector<double>> mypredict_probs=read_svm_output(myparavalue.pred_probs);
    
    cout<<"Get predict probs "<<endl;

    
  
    
    
    cout<<myparavalue.decoy_probs<<endl;

    vector<vector<double>> mydecoy_probs=read_svm_output(myparavalue.decoy_probs);
    
    
    cout<<"Get decoy probs"<<endl;
    
    
    cout<<"Get mypredict probs and decoy probs size "<<mypredict_probs.size()<<" "<<mydecoy_probs.size()<<endl;

    
    int can_col=myparavalue.can_col;
    
    int decoy_col=myparavalue.decoy_col;
    
    
    string decoy_tag=myparavalue.decoy_tag;
    
    
    
   // vector<int> mypredict_states=read_int(myparavalue.pred_states);
    
    
    
    
    
    //cout<<"Check reading size "<<mypredict_probs[0].size()<<" "<<mypredict_probs[1].size()<<endl;
    
    //cout<<mypredict_probs[0][0]<<" "<<mypredict_probs[1][0]<<" "<<mypredict_probs[2][0]<<endl;
    
    
    
    vector<string>myprotids=read_tsv(myparavalue.my_prots);
    
    //cout<<"get protids."<<endl;
    
    vector<int> mycan_protids=read_int(myparavalue.can_prots);
    
   // cout<<"get can_id"<<endl;
    
    vector<int> mydecoy_protids=read_int(myparavalue.decoy_prots);
    
    

    //cout<<"get the files"<<endl;
    
    
    vector<double>all_predict_probs(mypredict_probs.size());
    
    for (size_t i=0; i<mypredict_probs.size(); i++)
    {
        all_predict_probs[i]=mypredict_probs[i][can_col];
    }
    
    
    
    vector<double>all_decoy_probs(mydecoy_probs.size());
    
    for (size_t i=0; i<mydecoy_probs.size(); i++)
    {
        all_decoy_probs[i]=mydecoy_probs[i][decoy_col];
    }
    
    
    
    size_t numofprots=myprotids.size();
    
    vector<vector<double>> pred_probs_spec{myprotids.size()};
    
    for (size_t i=0; i<numofprots; i++)
    {
        for (size_t p=0; p<mycan_protids.size(); p++)
        {
            if (mycan_protids[p]==i)
            {
                pred_probs_spec[i].push_back(all_predict_probs[p]);
            }
        }
    }
    
    vector<vector<double>> decoy_probs_spec{myprotids.size()};
    
    for (size_t i=0; i<numofprots; i++)
    {
        
        for (size_t p=0; p<mydecoy_protids.size(); p++)
        {
            if (mydecoy_protids[p]==i)
            {
                decoy_probs_spec[i].push_back(all_decoy_probs[p]);
            }
        }
    }
    
    
    
    vector<vector<double>>each_prot_dfdr{999};
    
    for (int h=1; h<1000; h++)
    {
        double this_score=double(h)/double(1000);
        
        vector<double>at_each_dfdr(myprotids.size());
        
        for (size_t i=0; i<myprotids.size(); i++)
        {
            
            int numup=0;
            int numdown=0;
            
            for (size_t p=0; p<decoy_probs_spec[i].size(); p++)
            {
                if (decoy_probs_spec[i][p]>=this_score)
                {
                    numup=numup+1;
        
                }
                
            }
            
            
            for (size_t p=0; p<pred_probs_spec[i].size(); p++)
            {
                if (pred_probs_spec[i][p]>=this_score)
                {
                    numdown=numdown+1;
                }
                
            }
            
            if (numdown!=0 && decoy_probs_spec[i].size()!=0)
            {
                at_each_dfdr[i]=(double(numup)/double(numdown))*(double(pred_probs_spec[i].size())/double(decoy_probs_spec[i].size()));

            }
            else
            {
               at_each_dfdr[i]=-1;
            }
            
        
            
        }
        
        each_prot_dfdr[h-1]=at_each_dfdr;
        
    }
    
    
    ofstream outfile;
    
    cout<<"Give output "<<endl;
    
    outfile.open("prot_specific_der_" + decoy_tag + ".tsv");
    
    for (size_t i=0; i<each_prot_dfdr[0].size(); i++)
    {
        
        for (size_t p=0; p<999; p++)
        {
            outfile<<each_prot_dfdr[p][i]<<"\t";
        }
        
        outfile<<endl;
        
        
    }
    
    outfile.close();
    
    
    
    

    
    vector<double> all_dfdr(999);
    for (int h=1; h<1000; h++)
    {
        double this_score=double(h)/double(1000);
        
       
        int numup=0;
        int numdown=0;
        size_t all_can_size=0;
        size_t all_decoy_size=0;
        
        for (size_t i=0; i<myprotids.size(); i++)
        {
            
            for (size_t p=0; p<decoy_probs_spec[i].size(); p++)
            {
                if (decoy_probs_spec[i][p]>=this_score)
                {
                    numup=numup+1;
                    
                }
                
            }
            
            
            for (size_t p=0; p<pred_probs_spec[i].size(); p++)
            {
                if (pred_probs_spec[i][p]>=this_score)
                {
                    

                    numdown=numdown+1;
                }
                
            }
            
            all_can_size=all_can_size+pred_probs_spec[i].size();
            all_decoy_size=all_decoy_size+decoy_probs_spec[i].size();
            
        }
        
        
        
        
        if (numdown!=0)
        {
            
            //cout<<"monitor the ajustment size "<<all_can_size<<" "<<all_decoy_size<<endl;
            
            all_dfdr[h-1]=(double(numup)/double(numdown))*(double(all_can_size)/double(all_decoy_size));
            
        }
        else
        {
            all_dfdr[h-1]=-1;
        }
        
        
    }
    

    
    outfile.open("global_der_" + decoy_tag+ ".tsv");
    
    for (size_t i=0; i<all_dfdr.size(); i++)
    {
        outfile<<all_dfdr[i]<<endl;
    }
    
    outfile.close();
    
    
    cout<<"Done!"<<endl;
    
    
    
}

