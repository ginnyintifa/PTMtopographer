//
//  main.cpp
//  humansites
//
//  Created by Ginny Li on 15-7-28.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//



#include "site.hpp"
#include "aaprop.hpp"
#include "emisson.hpp" //not in use

#include "math.hpp"
#include "input_feature_generation.hpp"

int main(int argc, char *argv[])

{
    
    using namespace std;
    
    

    vector<string>myinput=read_input(argv[1]);
    vector<string>mypara=get_parameterlines(myinput);
    
    Para_values myparavalue=get_paravalues(mypara);
    vector<int> myaawhich=get_aawhich(myparavalue.pred_sites);
    
    string data_type=myparavalue.data_type;
    
    
    vector<prot_type> myproteins=read_fasta(myparavalue.fasta_name);//function using here is dependent on the format of file.
    

    cout<<"SEE how many proteins i have "<<myproteins.size()<<endl;
    

    
    vector<aa_properties> aminoacid_properties(24);
    
    assign_properties(aminoacid_properties);
    
    vector<double>each_hydro;
    vector<double>each_pka1;
    vector<double>each_helixpro;
    vector<double>each_stericpar;
    vector<double>each_polar;
    vector<double>each_volume;
    vector<double>each_sheetpro;
    

    assign_each_property(aminoacid_properties,each_hydro,each_pka1,each_helixpro,each_stericpar,each_polar,each_volume,each_sheetpro);

    
  
    vector<vector<string>>all_can_flank_windows{myproteins.size()};
    vector<vector<string>>all_decoy_flank_windows{myproteins.size()};
    vector<vector<string>>all_decoy_center_flank_windows{myproteins.size()};
    
    
    for (size_t i=0; i<myproteins.size(); i++)
    {
        vector<string>each_flank_window=get_windows_for_each(myproteins[i].second, myparavalue.flank_length, myparavalue.pred_sites);
        vector<string>each_decoy_window=get_decoy_for_each(myproteins[i].second, myparavalue.flank_length, myparavalue.pred_sites);
        
        vector<string> each_center_decoy_window=get_decoy_for_each(myproteins[i].second, 2, myparavalue.pred_sites); //get the center windows
        
        all_can_flank_windows[i]=each_flank_window;
        all_decoy_flank_windows[i]=each_decoy_window;
        all_decoy_center_flank_windows[i]=each_center_decoy_window;
        
    }
    

    
    
    vector<vector<vector<int>>>aa_compo_forall=get_allseqs__aavector(all_can_flank_windows, aminoacid_properties);
    
    vector<vector<vector<int>>>decoy_compo_forall=get_allseqs__aavector(all_decoy_flank_windows, aminoacid_properties);
    
    vector<vector<vector<int>>>decoy_center_compo_forall=get_allseqs__aavector(all_decoy_center_flank_windows, aminoacid_properties);
    
    
    for (size_t i=0; i<aa_compo_forall.size(); i++)
    {
         //cout<<i<<" "<<aa_compo_forall[i].size()<<endl;
    }
    
    
    vector<vector<vector<int>>>right_decoy_compo_forall{decoy_center_compo_forall.size()};
    vector<vector<vector<int>>>pure_decoy_compo_forall{decoy_center_compo_forall.size()};
    vector<vector<vector<int>>>near_decoy_compo_forall{decoy_center_compo_forall.size()};
    
    
    for (size_t i=0; i<decoy_center_compo_forall.size(); i++)
    {
        
        
        for (size_t p=0; p<decoy_center_compo_forall[i].size(); p++)
        {
            
            size_t insig=0;
            
            size_t pure_insig=0;
            
            for (size_t h=0; h<myaawhich.size(); h++)
            {
                if (decoy_center_compo_forall[i][p][myaawhich[h]]==0)
                {
                    insig=insig+1;
                    if (decoy_compo_forall[i][p][myaawhich[h]]==0)
                    {
                        pure_insig=pure_insig+1;
                    }
                    
                }
                
            }
            
            if (insig==myaawhich.size())
            {
                right_decoy_compo_forall[i].push_back(decoy_compo_forall[i][p]);
                
                
                if (pure_insig==myaawhich.size())
                {
                    pure_decoy_compo_forall[i].push_back(decoy_compo_forall[i][p]);
                }
                else
                {
                    near_decoy_compo_forall[i].push_back(decoy_compo_forall[i][p]);
                }
                
                
            }
            
            
            
        }
        
    }
    
    
    
    decoy_compo_forall=near_decoy_compo_forall; //Im using near decoy now
    
    
    
    
    
    vector<vector<int>>aa_vectorforall;//for the whole sequence
    
    for (size_t i=0; i<myproteins.size(); i++)
    {
        vector<int> aa_vectorforone;
        
        aa_vectorforone=grab_aa_vector(aminoacid_properties, myproteins[i].second);
        
        aa_vectorforall.push_back(aa_vectorforone);
    }
    
    
    
    vector<double>seqs_means_hydro=get_means(each_hydro,aa_vectorforall);
    
    vector<double>seqs_means_pka1=get_means(each_pka1,aa_vectorforall);
    
    vector<double>seqs_means_helixpro=get_means(each_helixpro,aa_vectorforall);
    
    vector<double>seqs_means_stericpar=get_means(each_stericpar,aa_vectorforall);
    
    vector<double>seqs_means_polar=get_means(each_polar,aa_vectorforall);
    
    vector<double>seqs_means_volume=get_means(each_volume,aa_vectorforall);
    
    vector<double>seqs_means_sheetpro=get_means(each_sheetpro,aa_vectorforall);
    
    
    
    //think about a way to have candidate/decoy/pure_decoy unified so I don't have to write so
    //many times
    //the key is to have the properties together
    
    
    
    vector<vector<double>> hydro_forall=get_allseqs_property(aa_compo_forall,each_hydro,seqs_means_hydro);//normalized
    
    vector<vector<double>> pka1_forall=get_allseqs_property(aa_compo_forall,each_pka1,seqs_means_pka1);
    
    vector<vector<double>> helixpro_forall=get_allseqs_property(aa_compo_forall,each_helixpro,seqs_means_helixpro);
    vector<vector<double>> stericpar_forall=get_allseqs_property(aa_compo_forall,each_stericpar,seqs_means_hydro);//normalized
    vector<vector<double>> polar_forall=get_allseqs_property(aa_compo_forall,each_polar,seqs_means_pka1);
    
    vector<vector<double>> volume_forall=get_allseqs_property(aa_compo_forall,each_volume,seqs_means_helixpro);
    vector<vector<double>> sheetpro_forall=get_allseqs_property(aa_compo_forall,each_sheetpro,seqs_means_helixpro);
    //////
    
    vector<vector<double>> decoy_hydro_forall=get_allseqs_property(decoy_compo_forall,each_hydro,seqs_means_hydro);//normalized
    
    vector<vector<double>> decoy_pka1_forall=get_allseqs_property(decoy_compo_forall,each_pka1,seqs_means_pka1);
    
    vector<vector<double>> decoy_helixpro_forall=get_allseqs_property(decoy_compo_forall,each_helixpro,seqs_means_helixpro);
    
    vector<vector<double>> decoy_stericpar_forall=get_allseqs_property(decoy_compo_forall,each_stericpar,seqs_means_hydro);//normalized
    
    vector<vector<double>> decoy_polar_forall=get_allseqs_property(decoy_compo_forall,each_polar,seqs_means_pka1);
    cout<<"see forall size "<<polar_forall.size()<<endl;
    
    vector<vector<double>> decoy_volume_forall=get_allseqs_property(decoy_compo_forall,each_volume,seqs_means_helixpro);
    
    vector<vector<double>> decoy_sheetpro_forall=get_allseqs_property(decoy_compo_forall,each_sheetpro,seqs_means_helixpro);

   //////
    vector<vector<double>> pure_decoy_hydro_forall=get_allseqs_property(pure_decoy_compo_forall,each_hydro,seqs_means_hydro);//normalized
    
    vector<vector<double>> pure_decoy_pka1_forall=get_allseqs_property(pure_decoy_compo_forall,each_pka1,seqs_means_pka1);
    
    vector<vector<double>> pure_decoy_helixpro_forall=get_allseqs_property(pure_decoy_compo_forall,each_helixpro,seqs_means_helixpro);
    
    vector<vector<double>> pure_decoy_stericpar_forall=get_allseqs_property(pure_decoy_compo_forall,each_stericpar,seqs_means_hydro);//normalized
    
    vector<vector<double>> pure_decoy_polar_forall=get_allseqs_property(pure_decoy_compo_forall,each_polar,seqs_means_pka1);
    cout<<"see forall size "<<polar_forall.size()<<endl;
    
    vector<vector<double>> pure_decoy_volume_forall=get_allseqs_property(pure_decoy_compo_forall,each_volume,seqs_means_helixpro);
    
    vector<vector<double>> pure_decoy_sheetpro_forall=get_allseqs_property(pure_decoy_compo_forall,each_sheetpro,seqs_means_helixpro);
  
    
    ofstream outfile;

    
//add the datatype variable to the name of the output file.
   
    
    outfile.open("can_sites_properties_"+ data_type +".tsv");
    
    
    outfile<<"Hydrophobicity"<<"\t"<<"Pka1"<<"\t"<<"Helix_frequency"<<"\t"<<"Steric_parameter"<<"\t"<<"Polar"<<"\t"<<"Volume"<<"\t"<<"Sheet_frequency"<<"\t";
    
    for (size_t k=0; k<20; k++)
    {
        outfile<<aminoacid_properties[k].aaname<<"\t";
    }
    
    
    
    outfile<<endl;
    
    
    
    for (size_t i=0 ; i<myproteins.size(); i++)
    {
        
        
        for (size_t p=0; p<aa_compo_forall[i].size(); p++)
        {
            outfile<<hydro_forall[i][p]<<"\t"<<pka1_forall[i][p]<<"\t"<<helixpro_forall[i][p]<<"\t"<<stericpar_forall[i][p]<<"\t"<<polar_forall[i][p]<<"\t"<<volume_forall[i][p]<<"\t"<<sheetpro_forall[i][p]<<"\t";
            
            for (size_t k=0; k<20; k++)
            {
                outfile<<aa_compo_forall[i][p][k]<<"\t";
            }
            
            outfile<<endl;
        }
    }
    
    outfile.close();
    
    
    
    
    
    outfile.open("decoy_sites_properties_"+ data_type +".tsv");
    
    
    outfile<<"Hydrophobicity"<<"\t"<<"Pka1"<<"\t"<<"Helix_frequency"<<"\t"<<"Steric_parameter"<<"\t"<<"Polar"<<"\t"<<"Volume"<<"\t"<<"Sheet_frequency"<<"\t";
    
    for (size_t k=0; k<20; k++)
    {
        outfile<<aminoacid_properties[k].aaname<<"\t";
    }
    
    
    
    outfile<<endl;
    
    
    
    for (size_t i=0 ; i<myproteins.size(); i++)
    {
        
        for (size_t p=0; p<decoy_compo_forall[i].size(); p++)
        {
            outfile<<decoy_hydro_forall[i][p]<<"\t"<<decoy_pka1_forall[i][p]<<"\t"<<decoy_helixpro_forall[i][p]<<"\t"<<decoy_stericpar_forall[i][p]<<"\t"<<decoy_polar_forall[i][p]<<"\t"<<decoy_volume_forall[i][p]<<"\t"<<decoy_sheetpro_forall[i][p]<<"\t";
            
            for (size_t k=0; k<20; k++)
            {
                outfile<<decoy_compo_forall[i][p][k]<<"\t";
            }
            
            outfile<<endl;
        }
    }
    
    outfile.close();
    
    
    outfile.open("pure_decoy_sites_properties_"+ data_type +".tsv");
    
    
    outfile<<"Hydrophobicity"<<"\t"<<"Pka1"<<"\t"<<"Helix_frequency"<<"\t"<<"Steric_parameter"<<"\t"<<"Polar"<<"\t"<<"Volume"<<"\t"<<"Sheet_frequency"<<"\t";
    
    for (size_t k=0; k<20; k++)
    {
        outfile<<aminoacid_properties[k].aaname<<"\t";
    }
    
    
    
    outfile<<endl;
    
    
    
    for (size_t i=0 ; i<myproteins.size(); i++)
    {
        
        for (size_t p=0; p<pure_decoy_compo_forall[i].size(); p++)
        {
            outfile<<pure_decoy_hydro_forall[i][p]<<"\t"<<pure_decoy_pka1_forall[i][p]<<"\t"<<pure_decoy_helixpro_forall[i][p]<<"\t"<<pure_decoy_stericpar_forall[i][p]<<"\t"<<pure_decoy_polar_forall[i][p]<<"\t"<<pure_decoy_volume_forall[i][p]<<"\t"<<pure_decoy_sheetpro_forall[i][p]<<"\t";
            
            for (size_t k=0; k<20; k++)
            {
                outfile<<pure_decoy_compo_forall[i][p][k]<<"\t";
            }
            
            outfile<<endl;
        }
    }
    
    outfile.close();

    
    
    

        cout<<myparavalue.backinfo<<endl;
        vector<string> mybackinfo=read_tsv(myparavalue.backinfo);
        
        cout<<mybackinfo.size()<<endl;
        for (size_t i=0; i<mybackinfo.size(); i++)
        {
            for (size_t p=0; p<mybackinfo[i].size(); p++)
            {
                int thesign=0;
                
                for (size_t q=0; q<myparavalue.pred_sites.size(); q++)
                {
                    if (mybackinfo[i][p]==myparavalue.pred_sites[q])
                    {
                        thesign=1;
                    }
                    
                    
                }
                
                if (thesign==0)
                {
                    mybackinfo[i][p]=toupper(mybackinfo[i][p]);
                }
            }
        }
        
        
        cout<<"get signs...\n";
        
        
        vector<string>mybackinfopeptides=mybackinfo;
        
        for (size_t i=0; i<mybackinfo.size(); i++)
        {
            for (size_t p=0; p<mybackinfo[i].size(); p++)
            {
                mybackinfopeptides[i][p]=toupper(mybackinfopeptides[i][p]);
            }
        }
        
        cout<<"get the peptides...\n";
        
    
        mybackinfopeptides=read_concise(mybackinfopeptides);
        
        
        vector<Pep_in_prot> mybacklinks=link_peptides_with_sites(mybackinfopeptides, mybackinfo,myparavalue.pred_sites);
        cout<<"get links...\n";
    
    
    
        pair<vector<prot_type>,vector<vector<Pep_in_prot>>>backppmap=get_mapps(myproteins, mybacklinks);
    

    
        cout<<"get maps...\n";
        vector<pair<prot_type,vector<size_t>>> mybackprotsitesall=get_all_seq_sites(backppmap);
        vector<pair<prot_type,vector<size_t>>> mybackneededprotsiteall=get_needed_myprot_sites_all(mybackprotsitesall);
    
    
    
        vector<pair<prot_type,vector<size_t>>> mybackuniquesites=get_unique_site(mybackneededprotsiteall);
    
    
    cout<<mybackuniquesites.size()<<endl;
    
    
    
    cout<<"see after psp mapping proteins "<<mybackuniquesites.size()<<endl;
    
    cout<<myparavalue.pred_sites[0]<<endl;
    
    vector<vector<int>>all_can_flank_states{myproteins.size()};
    vector<vector<size_t>>all_prot_can_pos{myproteins.size()};
    
    
    
    for (size_t i=0; i<myproteins.size(); i++)
    {
             for (size_t p=0; p<myproteins[i].second.size(); p++)
                {
                    for (size_t r=0; r<myparavalue.pred_sites.size(); r++)
                    {
                        if (myproteins[i].second[p]==toupper(myparavalue.pred_sites[r]))
                        {
                            all_can_flank_states[i].push_back(0);
                            
                            all_prot_can_pos[i].push_back(p);//have position recorded
                        }
                    }
                    
                    
                }
        
        
    }
    
    
    
    
  
    
    outfile.open("can_sites_position_"+ data_type +".tsv");
    
    for (size_t i=0; i<all_prot_can_pos.size(); i++)
    {
        for (size_t p=0; p<all_prot_can_pos[i].size(); p++)
        
        {
            outfile<<all_prot_can_pos[i][p]<<"\t";
        }
        
        outfile<<endl;
    }
    outfile.close();
    
  
    
    for (size_t i=0; i<myproteins.size(); i++)
    {
        for (size_t h=0 ; h<mybackuniquesites.size(); h++)
        {
            if (mybackuniquesites[h].first.first==myproteins[i].first)
            {
                
                for (size_t j=0; j<all_can_flank_states[i].size(); j++)
                {
                    
                    for (size_t p=0; p<mybackuniquesites[h].second.size(); p++)
                    {
                       
                        
                        if (all_prot_can_pos[i][j]==mybackuniquesites[h].second[p])
                        {
                            all_can_flank_states[i][j]=1;
                        }
                    }
                    
                    
                    
                }
                
                
            }
        }
        
      
    }
    
    
    outfile.open("head_annotation_"+ data_type +".tsv");
    
    
    outfile<<"protein"<<"\t"<<"protein_ID"<<"\t"<<"window"<<"\t"<<"position"<<"\t"<<"psp_state"<<"\t"<<"Hydrophobicity"<<"\t"<<"Pka1"<<"\t"<<"Helix_frequency"<<"\t"<<"Steric_parameter"<<"\t"<<"Polar"<<"\t"<<"Volume"<<"\t"<<"Sheet_frequency"<<"\t";
    
    for (size_t k=0; k<20; k++)
    {
        outfile<<aminoacid_properties[k].aaname<<"\t";
    }
    
    
    
    outfile<<endl;
    
    
    
    for (size_t i=0 ; i<myproteins.size(); i++)
    {
        
        //cout<<i<<" "<<myproteins[i].second.size()<<" aacompoforall_size "<<aa_compo_forall[i].size()<<endl;
        
        for (size_t p=0; p<aa_compo_forall[i].size(); p++)
        {
            outfile<<i+1<<"\t"<<myproteins[i].first<<"\t"<<all_can_flank_windows[i][p]<<"\t"<<all_prot_can_pos[i][p]+1<<"\t"<<all_can_flank_states[i][p]<<"\t"<<hydro_forall[i][p]<<"\t"<<pka1_forall[i][p]<<"\t"<<helixpro_forall[i][p]<<"\t"<<stericpar_forall[i][p]<<"\t"<<polar_forall[i][p]<<"\t"<<volume_forall[i][p]<<"\t"<<sheetpro_forall[i][p]<<"\t";
            
            for (size_t k=0; k<20; k++)
            {
                outfile<<aa_compo_forall[i][p][k]<<"\t";
            }
            
            outfile<<endl;
        }
    }
    
    outfile.close();
    

    
    
    outfile.open("protein_head_annotation_"+ data_type +".tsv");
    
    outfile<<"protein"<<"\t"<<"protein_ID"<<"\t"<<"protein_length"<<"\t"<<"#candidate_windows"<<"\t"<<"#decoy_windows"<<"\t"<<"#pure_decoy_windows"<<"\t"<<"ave_Hydrophobicity"<<"\t"<<"ave_Pka1"<<"\t"<<"ave_Helix_frequency"<<"\t"<<"ave_Steric_parameter"<<"\t"<<"ave_Polar"<<"\t"<<"ave_Volume"<<"\t"<<"ave_Sheet_frequency"<<"\t";
    
    for (size_t k=0; k<20; k++)
    {
        outfile<<aminoacid_properties[k].aaname<<"\t";
    }
    
    outfile<<endl;

    
    for (size_t i=0; i<myproteins.size(); i++)
    {
        outfile<<i+1<<"\t"<<myproteins[i].first<<"\t"<<myproteins[i].second.size()<<"\t"<<aa_compo_forall[i].size()<<"\t"<<decoy_compo_forall[i].size()<<"\t"<<pure_decoy_compo_forall[i].size()<<"\t"<<seqs_means_hydro[i]<<"\t"<<seqs_means_pka1[i]<<"\t"<<seqs_means_helixpro[i]<<"\t"<<seqs_means_stericpar[i]<<"\t"<<seqs_means_polar[i]<<"\t"<<seqs_means_volume[i]<<"\t"<<seqs_means_sheetpro[i]<<"\t";
        
        
        for (size_t k=0; k<20; k++)
        {
            outfile<<aa_vectorforall[i][k]<<"\t";
        }
        
        outfile<<endl;
        
    }
    
    outfile.close();
    
    
    
    
    
    
    
    
    outfile.open("can_sites_states_"+ data_type +".tsv");
    
    for (size_t i=0; i<all_can_flank_states.size(); i++)
    {
        for (size_t p=0; p<all_can_flank_states[i].size(); p++)
        {
            outfile<<all_can_flank_states[i][p]<<endl;
        }
    }
    
    outfile.close();
    
    
    outfile.open("my_proteins_"+ data_type +".tsv");
    for (size_t i=0; i<myproteins.size(); i++)
    {
        outfile<<myproteins[i].first<<endl;
    }
    outfile.close();
    
    
    outfile.open("can_sites_prots_"+ data_type +".tsv");
    
    
    for (size_t i=0; i<myproteins.size(); i++)
    {
        for (size_t p=0; p<all_can_flank_states[i].size(); p++)
        {
            outfile<<i<<endl;
        }
        
       
    }
    
    outfile.close();
    
    outfile.open("decoy_sites_prots_"+ data_type +".tsv");
    for (size_t i=0; i<myproteins.size(); i++)
    {
        for (size_t p=0; p<decoy_compo_forall[i].size(); p++)
        {
            outfile<<i<<endl;
        }
    }
    
    
    outfile.close();
    
    
    outfile.open("pure_decoy_sites_prots_"+ data_type +".tsv");
    for (size_t i=0; i<myproteins.size(); i++)
    {
        for (size_t p=0; p<pure_decoy_compo_forall[i].size(); p++)
        {
            outfile<<i<<endl;
        }
    }
    
    
    outfile.close();

    
    
    
    
    
    outfile.open("can_svm_input_"+ data_type +".tsv");
    
    
    for (size_t i=0 ; i<myproteins.size(); i++)
    {
        
        for (size_t p=0; p<aa_compo_forall[i].size(); p++)
        {
            outfile<<all_can_flank_states[i][p]<<"\t"<<"1:"<<hydro_forall[i][p]<<"\t"<<"2:"<<pka1_forall[i][p]<<"\t"<<"3:"<<helixpro_forall[i][p]<<"\t"<<"4:"<<stericpar_forall[i][p]<<"\t"<<"5:"<<polar_forall[i][p]<<"\t"<<"6:"<<volume_forall[i][p]<<"\t"<<"7:"<<sheetpro_forall[i][p]<<"\t";
            
            for (size_t k=0; k<20; k++)
            {
                outfile<<k+8<<":"<<aa_compo_forall[i][p][k]<<"\t";
            }
            
            outfile<<endl;
        }
    }
    
    outfile.close();
    
    
    
    
    
    outfile.open("decoy_svm_input_"+ data_type +".tsv");
    
   
    for (size_t i=0 ; i<myproteins.size(); i++)
    {
        
        
        for (size_t p=0; p<decoy_compo_forall[i].size(); p++)
        {
            outfile<<0<<"\t"<<"1:"<<decoy_hydro_forall[i][p]<<"\t"<<"2:"<<decoy_pka1_forall[i][p]<<"\t"<<"3:"<<decoy_helixpro_forall[i][p]<<"\t"<<"4:"<<decoy_stericpar_forall[i][p]<<"\t"<<"5:"<<decoy_polar_forall[i][p]<<"\t"<<"6:"<<decoy_volume_forall[i][p]<<"\t"<<"7:"<<decoy_sheetpro_forall[i][p]<<"\t";
            
            for (size_t k=0; k<20; k++)
            {
                outfile<<k+8<<":"<<decoy_compo_forall[i][p][k]<<"\t";
            }
            
            outfile<<endl;
        }
    }
    
    outfile.close();
    
    
    outfile.open("pure_decoy_svm_input_"+ data_type +".tsv");
    
    
    for (size_t i=0 ; i<myproteins.size(); i++)
    {
        
       
        
        for (size_t p=0; p<pure_decoy_compo_forall[i].size(); p++)
        {
            outfile<<0<<"\t"<<"1:"<<pure_decoy_hydro_forall[i][p]<<"\t"<<"2:"<<pure_decoy_pka1_forall[i][p]<<"\t"<<"3:"<<pure_decoy_helixpro_forall[i][p]<<"\t"<<"4:"<<pure_decoy_stericpar_forall[i][p]<<"\t"<<"5:"<<pure_decoy_polar_forall[i][p]<<"\t"<<"6:"<<pure_decoy_volume_forall[i][p]<<"\t"<<"7:"<<pure_decoy_sheetpro_forall[i][p]<<"\t";
            
            for (size_t k=0; k<20; k++)
            {
                outfile<<k+8<<":"<<pure_decoy_compo_forall[i][p][k]<<"\t";
            }
            
            outfile<<endl;
        }
    }
    
    outfile.close();
    

    
    
    
    
    
    
    
    
}

