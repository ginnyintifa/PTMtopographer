//
//  site.hpp
//  humansites
//
//  Created by Ginny Li on 15-7-29.
//  Copyright (c) 2015年 Ginny. All rights reserved.
//





#ifndef site_hpp
#define site_hpp


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
#include <sstream>

//http://www.boost.org/doc/libs/master/libs/algorithm/doc/html/algorithm/Searching.html
// #include "boyer_moore_horspool.hpp"
#include <boost/algorithm/searching/boyer_moore_horspool.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>


template <class CharT, class Traits>
void check_init_istream(::std::basic_istream<CharT, Traits>& is)
{
	is.exceptions(::std::ios_base::failbit);
	is.exceptions(::std::ios_base::badbit);
}





typedef std::pair<std::string,std::string> prot_type;





/*
struct Site_str{
    std::vector<std::pair<char,size_t>> site_position;
    //double logfoldchange;
};
*/




struct Pep_in_prot{
    std::string pep_seq;
    std::size_t mpos;//the position of the peptide in a mapping ;
    std::vector<size_t> sites_positions;
    
};




struct search_res_type{
	std::size_t prot_idx; // protein index
	std::size_t pep_pos;// position where peptide is found
};



std::vector<prot_type>
read_fasta(const std::string& fn){
	std::ifstream is(fn); check_init_istream(is);
	std::vector<prot_type> ret;
    std::string header, seq;
	is.ignore(1);// omit the first ">" character
	for(;;){
		getline(is,header);// read the header line
		getline(is,seq,'>');// read the lines representing a sequence
		if(!is) break;
		const auto endit=std::remove(seq.begin(),seq.end(),'\n');// remove all newlines from the sequence
		ret.emplace_back(
                         header,// ID
                         std::string(seq.begin(),endit)// sequence
                         );
	}
	return ret;
}


std::vector<prot_type>
read_non_isoform(const std::string& fn)
{
    using namespace std;
    ifstream myfile(fn); check_init_istream(myfile);
    
    vector<prot_type> prots;
    string header,seq;
    for(;;)
    {
        getline(myfile,header);
        getline(myfile,seq);
        
        if(!myfile) break;
        
        prots.emplace_back(header, string(seq.begin(),seq.end()));
    }
    
    return prots;
    
}




std::vector<int> read_int(const std::string& fn)
{
    using namespace std;
    ifstream ifs(fn);
    check_init_istream(ifs);
    vector<int> v;
    string line;
    while(getline(ifs,line))
        v.push_back(
                    line.empty()?
                    numeric_limits<int>::quiet_NaN():
                    stod(line));
    return v;
}








std::vector<std::string> read_tsv(const std::string& fn)
{
	using namespace std;
	ifstream myfile(fn); check_init_istream(myfile);
    
	vector<std::string> cpeptides;
    std::string line;
	while(getline (myfile,line))
		cpeptides.emplace_back(line);
    
	return cpeptides;
}



//put the ratio missing value as NAN.
//not being used anymore
std::vector<double> read_tsv_double(const std::string& fn)
{
	using namespace std;
	ifstream ifs(fn);
    check_init_istream(ifs);
	vector<double> v;
    string line;
	while(getline(ifs,line))
		v.push_back(
                    line.empty()?
                    numeric_limits<double>::quiet_NaN():
                    stod(line));
	return v;
}



std::vector<std::vector<double>> read_svm_output(const std::string& fn)
{
    using namespace std;
    
    cout<<"first check "<<endl;
    
    ifstream ifs(fn);
    
    cout<<"pass first"<<endl;
    check_init_istream(ifs);
    
    cout<<"start reading "<<endl;
    
    string dummyLine;
    getline(ifs,dummyLine);
    
    cout<<dummyLine<<endl;
    
    vector<vector<double>> v;
    
    string line;
    while (getline(ifs,line))
    {
        istringstream is( line );
        v.push_back(vector<double>(istream_iterator<double>(is), istream_iterator<double>()));
    }

    return v;
    
}






std::vector<std::vector<int>> read_tsv_int(const std::string& fn)
{
    using namespace std;
    ifstream ifs(fn);
    check_init_istream(ifs);
    vector<vector<int>> v;
    string line;
    while (getline(ifs,line))
    {
        istringstream is( line );
        v.push_back(vector<int>(istream_iterator<int>(is), istream_iterator<int>()));
    }
    
    return v;
    
}


std::vector<std::vector<double>> read_vector_double(const std::string& fn)
{
    using namespace std;
    ifstream ifs(fn);
    check_init_istream(ifs);
    vector<vector<double>> v;
    string line;
    while (getline(ifs,line))
    {
        istringstream is( line );
        v.push_back(vector<double>(istream_iterator<double>(is), istream_iterator<double>()));
    }
    
    return v;
    
}





//a function that read the emission into structure wanted

std::vector<std::vector<std::pair<double,double>>>get_emi_structure(const std::vector<double>&readtsv)
{
    using namespace std;
    
    size_t lenofeach=readtsv.size()/4;
    
    vector<vector<pair<double,double>>>prop_emission{2};
    
    for (size_t i=0 ; i<lenofeach; i++)
    {
        pair<double,double>state0_emis;
        state0_emis.first=readtsv[i];
        state0_emis.second=readtsv[i+2*lenofeach];
        
        pair<double,double>state1_emis;
        state1_emis.first=readtsv[i+lenofeach];
        state1_emis.second=readtsv[i+3*lenofeach];
        
        prop_emission[0].push_back(state0_emis);
        prop_emission[1].push_back(state1_emis);
    }
    
    return prop_emission;
    
}

std::vector<std::vector<double>>get_aafreq_structure(const std::vector<double>&readtsv)
{
    using namespace std;
    size_t lenofeach=readtsv.size()/2;
   // cout<<lenofeach<<endl;
    
    vector<vector<double>>aa_freq{2};
    
    for (size_t i=0; i<lenofeach; i++)
    {
        aa_freq[0].push_back(readtsv[i]);
        aa_freq[1].push_back(readtsv[i+lenofeach]);
    }
    
    return aa_freq;
    
}


/////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////get concised proteins and peptides///////////////////////////////


std::vector<prot_type>read_concise_pair(const std::vector<prot_type>& things)
{
	using namespace std;
	unordered_set<prot_type,boost::hash<prot_type>> s(things.begin(), things.end());
	return vector<prot_type> (s.begin(), s.end());
}


std::vector<size_t>read_concise_position(const std::vector<size_t>& positions)
{
    using namespace std;
    unordered_set<size_t, boost::hash<size_t>>  s(positions.begin(),positions.end());
    return vector<size_t>(s.begin(),s.end());
}





std::vector<std::string>read_concise(const std::vector<std::string>& mystrings)
{
    using namespace std;
    unordered_set<string, boost::hash<string>> s(mystrings.begin(),mystrings.end());
    return vector<string> (s.begin(),s.end());
}




std::vector<size_t>read_concise_site(const std::vector<size_t>& oneprotsites)
{
    using namespace std;
    
    unordered_set<size_t,boost::hash<size_t>> s(oneprotsites.begin(),oneprotsites.end());
    
    return vector<size_t> (s.begin(),s.end());
}



//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////divide_into windows///////////////////////////////////////

//to be modified


//only for candidate sites

std::vector<std::string> get_windows_for_each(std::string my_sequence,int flank_length, std::string pred_sites)
{
    using namespace std;
    
    vector<string> my_windows;
    
    
    for (size_t i=0; i<my_sequence.size(); i++)
    {
       // cout<<"start to process this amino acid "<<i<<" ";
        
        for (size_t h=0; h<pred_sites.size(); h++)
        {
            if (my_sequence[i]==toupper(pred_sites[h]))
            {
         //       cout<<"enter selection "<<i<<" ";
                
                string one_window;
                size_t start_pos;
                size_t end_pos;
                
                start_pos=(((int(i)-flank_length)>0)?(int(i)-flank_length):0);
           //     cout<<"start at "<<start_pos<<" ";
                end_pos=(((i+flank_length)>(my_sequence.size()-1))?(my_sequence.size()-1):(i+flank_length));
             //   cout<<"end at "<<end_pos<<" ";
                one_window=my_sequence.substr(start_pos,(end_pos-start_pos)+1);
                
                my_windows.push_back(one_window);
                
            }
            
            //cout<<i<<" size "<<my_windows.size()<<endl;
        }
    }

    
    return my_windows;
}





// I need to change the way of selecting decoy
std::vector<std::string> get_decoy_for_each(std::string my_sequence,int flank_length, std::string pred_sites)
{
    using namespace std;
    
    vector<string> my_decoy;
    
    
    for (size_t i=0; i<my_sequence.size(); i++)
    {
        // cout<<"start to process this amino acid "<<i<<" ";
        
        for (size_t h=0; h<pred_sites.size(); h++)
        {
            if (my_sequence[i]!=toupper(pred_sites[h]))
            {
                //       cout<<"enter selection "<<i<<" ";
                
                string one_window;
                size_t start_pos;
                size_t end_pos;
                
                start_pos=(((int(i)-flank_length)>0)?(int(i)-flank_length):0);
                //     cout<<"start at "<<start_pos<<" ";
                end_pos=(((i+flank_length)>(my_sequence.size()-1))?(my_sequence.size()-1):(i+flank_length));
                //   cout<<"end at "<<end_pos<<" ";
                one_window=my_sequence.substr(start_pos,(end_pos-start_pos)+1);
                
                my_decoy.push_back(one_window);
                
            }
            
            //cout<<i<<" size "<<my_windows.size()<<endl;
        }
    }
    
    
    return my_decoy;
}




std::vector<std::vector<std::string>>get_windows_for_all(std::vector<std::pair<prot_type,std::vector<size_t>>> mysinglesites,int flank_length,std::string pred_sites)
{
    using namespace std;
    
    vector<vector<string>> window_forall;
    
    for (size_t i=0; i<mysinglesites.size(); i++)
        
    {
       // cout<<"protein "<<i<<" "<<mysinglesites[i].first.second<<endl;
        
        vector<string> windows=get_windows_for_each(mysinglesites[i].first.second,flank_length,pred_sites);
        window_forall.push_back(windows);
    }
    
    return window_forall;
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////link with sites/////////////////////////////////////////////

//do the link between sites and peptides

//ratio should never be there now

//I need to rewrite the whole thing yes



std::vector<Pep_in_prot> link_peptides_with_sites(const std::vector<std::string>& myconcised_peptides, std::vector<std::string> mypeptides_sites,std::string pred_sites)
{
    using namespace std;
    
    vector<Pep_in_prot>mylinks(myconcised_peptides.size());
    
    
    //vector<pair<string,Site_str>> all_sites_positions;
    
    vector<pair<string, vector<size_t>>> all_pep_sites;
    
    for (size_t p=0;p<mypeptides_sites.size();p++)
    {
       
        pair<string, vector<size_t>>one_pep_sites;
        
        
        for (size_t m=0; m<mypeptides_sites[p].size(); m++)
        {
            
            //cout<<"target "<<pred_sites<<endl;
            
            for (size_t f=0; f<pred_sites.size(); f++)
            {
                

                
                if (mypeptides_sites[p][m]==pred_sites[f])//I get a problem here
                {
                    size_t asite;
                    
                    //cout<<"get an incidence "<<p<<" "<<m<<endl;
                    
                    asite=m;
                 
                    one_pep_sites.second.push_back(asite);
                    //I just move the toupper outside of if
                }
                

            }
            
            
            mypeptides_sites[p][m]=toupper(mypeptides_sites[p][m]);
            

        }
        
        
        one_pep_sites.first=mypeptides_sites[p];
        
        all_pep_sites.push_back(one_pep_sites);
        
    }
    
    
    
    for (size_t i=0; i<myconcised_peptides.size(); i++)
    {
        for (size_t p=0; p<all_pep_sites.size(); p++)
        {
            if (all_pep_sites[p].first==myconcised_peptides[i])
            {
                for (size_t h=0; h<all_pep_sites[p].second.size(); h++)
                {
                    mylinks[i].sites_positions.push_back(all_pep_sites[p].second[h]);
                }
                
            }
       
        }
        
        mylinks[i].pep_seq=myconcised_peptides[i];
        mylinks[i].mpos=0;

    }
    
    
    return mylinks; //now you mayfind repeated site_type in one pep_in_prot(concise it later)
}



//.mpos is the place I put position

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
std::vector<search_res_type>
query_pep(
          const std::string& ncpeptide,
          const std::vector<prot_type>& newproteins
          )
{
	using namespace std;
	vector<search_res_type> index_pos;
		const auto search_func=boost::algorithm::make_boyer_moore_horspool(ncpeptide);
	for(size_t i=0;i<newproteins.size();++i) {
		const auto& prot_seq=newproteins[i].second;
        for(auto it=prot_seq.begin();;++it){
			it=search_func(it,prot_seq.end());
			if(it==prot_seq.end())// if not found
				break;
			search_res_type res;
			res.prot_idx=i;
			res.pep_pos=it-prot_seq.begin();
			index_pos.push_back(res);
		}
	}
	return index_pos;
}



std::pair<std::vector<prot_type>,std::vector<std::vector<Pep_in_prot>>>
get_mapps(const std::vector<prot_type>& myproteins,const std::vector<Pep_in_prot>& mylinks)
{
	using namespace std;

    
	vector<vector<search_res_type>> all_index_pos{mylinks.size()};
    
	cerr<<"searching...\n";
    
  #pragma omp parallel for
    
	for(size_t q=0;q<mylinks.size();q++)
    {
		all_index_pos[q]= query_pep(mylinks[q].pep_seq, myproteins);
    }
    
    
	cerr<<"to get all the ids and protein seqs which have peptides mapped on...\n";
    
    
	vector<vector<Pep_in_prot>> peps_sites_positions(myproteins.size());
    
	for (size_t r=0;r<all_index_pos.size();r++)
		BOOST_FOREACH(const auto& search, all_index_pos[r]){
			Pep_in_prot pep;
            
            
			pep.pep_seq=mylinks[r].pep_seq;
            pep.mpos=search.pep_pos;
            pep.sites_positions=mylinks[r].sites_positions;
			
			peps_sites_positions[search.prot_idx].push_back(pep);
		}
    
    cout<<"when doing the mapps "<<peps_sites_positions.size()<<endl;
	///// filter out proteins without peptide match
	vector<prot_type> myproteins2;
	vector<vector<Pep_in_prot>> peps_sites_positions2;
    
	for (size_t r=0;r<myproteins.size();r++)
		if(!peps_sites_positions[r].empty())//if not empty
        {
			myproteins2.push_back(myproteins[r]);
			peps_sites_positions2.push_back(peps_sites_positions[r]);
		}
    
  	return make_pair(myproteins2, peps_sites_positions2);
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::pair<prot_type,std::vector<size_t>>>get_all_seq_sites(const std::pair<std::vector<prot_type>,std::vector<std::vector<Pep_in_prot>>>&mymapps)
{
    using namespace std;
    vector<pair<prot_type,vector<size_t>>> myprot_sites_all;//for all the sequences
    
    for(size_t i=0;i<mymapps.first.size();i++)
    {
        pair<prot_type,vector<size_t>> myprot_sites_eachprot;
        
        myprot_sites_eachprot.first=mymapps.first[i];
        
        
        
        
        for (size_t p=0; p<mymapps.second[i].size(); p++)//the pth mapping of protein i;
        {
            size_t baseposition=mymapps.second[i][p].mpos;
            
            for(size_t k=0; k<mymapps.second[i][p].sites_positions.size();k++)
            {
                //double thelogfoldchange=mymapps.second[i]
                
                //[p].sites_positions_logfoldchange[k].logfoldchange;
               
                //something very wrong here!!!
                
                //about site_position
                
                
                //for(size_t h=0;h<mymapps.second[i][p].sites_positions[k].site_position.size();h++)
                //{
                    size_t myonesite;
                    //myonesite.aa=mymapps.second[i][p].sites_positions[k].aa;
                    myonesite=mymapps.second[i][p].sites_positions[k]+baseposition;
                    //myonesite.logfoldchange=thelogfoldchange;
                    
                    myprot_sites_eachprot.second.push_back(myonesite);
                    
                }
                
            }
            
            
        
        
        myprot_sites_all.push_back(myprot_sites_eachprot);
    }
    return myprot_sites_all;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool sortstruct(const size_t & asite, const size_t & bsite)
{
    return asite < bsite;
}


bool sortpair(const std::pair<prot_type, std::vector<size_t>> & asinglesites, const std::pair<prot_type, std::vector<size_t>> & bsinglesites)
{
    return asinglesites.second.size()>bsinglesites.second.size();
}


bool sortdoublepair(const std::pair<double, double> &aemi,const std::pair<double,double>& bemi )
{
    return aemi.first<bemi.first;
}



//I need to abandon the protein sequences which have no sites on




std::vector<std::pair<prot_type,std::vector<size_t>>>get_needed_myprot_sites_all(std::vector<std::pair<prot_type,std::vector<size_t>>> myprot_sites_all)
{
    using namespace std;
    
    vector<pair<prot_type,vector<size_t>>>myneeded_protsites_all;
    
    for(size_t i=0;i<myprot_sites_all.size();i++)
    {
        if(myprot_sites_all[i].second.size()>0)
            myneeded_protsites_all.push_back(myprot_sites_all[i]);
            }
    
    return myneeded_protsites_all;
    
    
    
}

//now the problem is my windows with sites only appear in the first half

//I should debug now



std::vector<std::pair<prot_type,std::vector<size_t>>>get_unique_site(std::vector<std::pair<prot_type,std::vector<size_t>>>myprot_sites_all)
{
    using namespace std;
    vector<pair<prot_type,vector<size_t>>> all_unique_site(myprot_sites_all.size());
    
    for (size_t i=0; i<myprot_sites_all.size(); i++)
    {
        //because the function is down now!!!
     
        vector<size_t>concise_protsites=myprot_sites_all[i].second;
        
        all_unique_site[i].first=myprot_sites_all[i].first;
        all_unique_site[i].second=read_concise_site(concise_protsites);
        
        
        for (size_t p=0; p<all_unique_site[i].second.size(); p++)
        {
           // cout<<" protein "<<i<<" site "<<p<<" which "<<all_unique_site[i].second[p]<<endl;
        }
        
        
        
    }
    
    return all_unique_site;
    
}



std::vector<std::pair<prot_type,std::vector<std::vector<size_t>>>>get_unique_site_inblock(const std::vector<std::pair<prot_type,std::vector<size_t>>>& all_unique_site,int flank_length, std::string pred_sites)
{
    using namespace std;
    
    vector<pair<prot_type,vector<vector<size_t>>>> all_uniquesites_inblock(all_unique_site.size());//the size of protein sequence.
    
    for (size_t i=0;i<all_unique_site.size();i++)
    {
        all_uniquesites_inblock[i].first=all_unique_site[i].first;
        
        size_t seqlen=all_unique_site[i].first.second.size();
        
        
        vector<pair<size_t,size_t>>window_start_end;
        
        
        for (size_t j=0; j<seqlen; j++)
        {
            for (size_t h=0; h<pred_sites.size(); h++)
            {
                if (all_unique_site[i].first.second[j]==toupper(pred_sites[h]))
                {
                    pair<size_t,size_t> one_start_end;
                    
                    one_start_end.first=(((int(j)-flank_length)>0)?(int(j)-flank_length):0);
                    one_start_end.second=(((j+flank_length)>(seqlen-1))?(seqlen-1):(j+flank_length));
                    
                    window_start_end.push_back(one_start_end);
                    
                }
            }
        }

        
        size_t numblock=window_start_end.size();
        
    
        vector<vector<size_t>>blocks_in_seq{numblock};
        

        
        for(size_t p=0;p<all_unique_site[i].second.size();p++)
            
        {
            size_t pos=all_unique_site[i].second[p];
            
            //cout<<"protein "<<i<<" site "<<p<<" position "<<pos<<endl;
            
            for (size_t h=0; h<numblock; h++)
            {
               // cout<<"block "<<h<<" start "<<window_start_end[h].first<<" end "<<window_start_end[h].second<<endl;
                
                if ((pos>=window_start_end[h].first)&(pos<=window_start_end[h].second))
                {
                   // cout<<"enter selection "<<h<<" position "<<pos<<endl;
                    
                    blocks_in_seq[h].push_back(all_unique_site[i].second[p]);
                }
                
                //cout<<"************************************protein "<<i<<" block "<<h<<" size "<<blocks_in_seq[h].size()<<endl;

            }
            
        }
        
        
        all_uniquesites_inblock[i].second=blocks_in_seq;
    }

    
            return all_uniquesites_inblock;
}




//have a look

std::vector<std::pair<prot_type, std::vector<int>>>get_blocks_withsites_all(std::vector<std::pair<prot_type,std::vector<std::vector<size_t>>>>all_uniquesites_inblock)
{
    using namespace std;
    
    vector<pair<prot_type,vector<int>>>all_blocks_withsites(all_uniquesites_inblock.size());
    
    for (size_t i=0; i<all_uniquesites_inblock.size(); i++)
    {
        all_blocks_withsites[i].first=all_uniquesites_inblock[i].first;
        
        vector<int> sitesblock;
        
        for (size_t p=0; p<all_uniquesites_inblock[i].second.size(); p++)
        {
            
            
            if(all_uniquesites_inblock[i].second[p].size()>0)
            {
                
                
                sitesblock.push_back(p);
            }
            
        }
        
        all_blocks_withsites[i].second=sitesblock;
    }
    
    return all_blocks_withsites;
        
    
    
    
}





//now I will make something that can be used to analyze the state of each window.


std::vector<std::vector<int>>get_states(std::vector<std::pair<prot_type,std::vector<std::vector<size_t>>>> sites_inwindows,std::string pred_sites)
{
    using namespace std;
    
    vector<vector<int>>windows_states;
    
    //cout<<"something wrong here"<<endl;
    
    //cout<<"the size "<<sites_inwindows.size()<<endl;
    
    for (size_t i=0; i<sites_inwindows.size(); i++)
    {
        vector<int>windows_states_seq(sites_inwindows[i].second.size(),0);
        //cout<<i<<endl;
        
        for (size_t p=0; p<sites_inwindows[i].second.size(); p++)
        {
            
            //cout<<"wrong "<<p<<endl;
            if (sites_inwindows[i].second[p].size()==0)
            {
                windows_states_seq[p]=0;
            }
            
            else
            {
                
                
                windows_states_seq[p]=1;
                
                /*
                for (size_t h=0; h<sites_inwindows[i].second[p].size(); h++)
                {
                 
                    for (size_t f=0; f<pred_sites.size(); f++)
                    {
                        if (sites_inwindows[i].first.second[sites_inwindows[i].second[p][h]]==toupper(pred_sites[f]))//'k')//input parameter come may also be problematic
                        {
                            windows_states_seq[p]=1;
                        }

                    }
                    
                }
                */
                
            }
            
          //  cout<<i<<" "<<p<<" $$$$$ "<<windows_states_seq[p]<<endl;
    
        }
        
        
        windows_states.push_back(windows_states_seq);
        
        
        
        
    }
    
  //  cout<<"finish here...\n";
    
    return windows_states;

}


//NOTE: for some protein sequence, for a certain state, there may not be any window in this particular state, so , there will be some 0s and nans when dividing
//NOTE:this is for non overlapping cases only, if there are overlapps I have to add another case.


std::vector<std::vector<std::string>>marked_windows(const std::vector<std::pair<prot_type,std::vector<std::vector<size_t>>>> &window_singlesites_forall,std::vector<std::vector<std::string>>window_forall,int flank_length,std::string pred_sites)
{
    using namespace std;
    
    vector<vector<string>>marked_window_forall;
    marked_window_forall=window_forall;
    
    

  
        for (size_t i=0; i<window_singlesites_forall.size();i++)
        {
            
            vector<pair<size_t,size_t>>window_start_end;
            
            size_t seqlen=window_singlesites_forall[i].first.second.size();
            
            
            for (size_t j=0; j<seqlen; j++)
            {
                for (size_t h=0; h<pred_sites.size(); h++)
                {
                    if (window_singlesites_forall[i].first.second[j]==toupper(pred_sites[h]))
                    {
                        pair<size_t,size_t> one_start_end;
                        
                        one_start_end.first=(((int(j)-flank_length)>0)?(int(j)-flank_length):0);
                        one_start_end.second=(((j+flank_length)>(seqlen-1))?(seqlen-1):(j+flank_length));
                        
                        window_start_end.push_back(one_start_end);
                        
                    }
                }
            }
            
            
            
            
            for (size_t p=0; p<window_singlesites_forall[i].second.size(); p++)//for each block
            {
               // cout<<i<<" how many sites in each block "<<p<<" "<<window_singlesites_forall[i].second[p].size()<<endl;
                if (window_singlesites_forall[i].second[p].size()>0)
                {
                    
                    for (size_t h=0; h<window_singlesites_forall[i].second[p].size(); h++)//for each site in this block
                    {
                        
                        size_t pos=window_singlesites_forall[i].second[p][h];
                       
                        
                        pos=pos-window_start_end[p].first;
                        
                        marked_window_forall[i][p][pos]=tolower(window_forall[i][p][pos]);
                        
                        // cout<<p<<" position of each site "<<window_singlesites_forall[i].second[p][h].position<<endl;
                    }
                }
            }
        }

    
    
    return marked_window_forall;
    
    
    
    
    
}







#endif
