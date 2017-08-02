# PTMtopographer user manual

Ginny ginnyli056@gmail.com / A0123847@u.nus.edu

This manual explains the workflow to run PTMtopographer, with an example of ubiquitination analysis(two-fold cross validation) on a set of target proteins. In this tutorial, we will use the Ubi_K-peptides from the PhosphoSitePlus database to train a random forest classifier, and compute prediction scores (the probability of having a Ubi_K) for candidate windows and decoy windows. The tutorial generates site-level annotation as well as protein-level annotation output.

Note that the input parameter files we provide is ready to use if you follow the tutorial.




## Step 0. Installation
For users with Mac OS X, GNU gcc compilers must be updated according to their OS versions: see http://hpc.sourceforge.net/. This can be done by

```
gunzip gcc-x-bin.tar.gz
sudo tar -xvf gcc-x-bin.tar -C/
```
Following this, Makefile is provided for PTMtopographer:

```
make
```
should compile the source and create three binary files in the /PTMtopographer/bin/.

After compiling, direct to bin for running：

```
cd bin
```



## Step 1. Preparing input files
* Fasta file containing training and test protein sequences (e.g. Uniprot fasta file).
* Training peptide sequences (e.g. PhosphoSitePlus) with PTM sites marked in lower case letters 

### 1A. Fasta file of sequences. 
One is training sequences, the other is test sequences. For example:

```
>tr|A0A024QZ18|A0A024QZ18_HUMAN PDZ domain containing, X chromosome, isoform CRA_a ...
MPLLWITGPRYHLILLSEASCLRANYVHLCPLFQHRWLETCNAPPQLIQGKARSAPKPSQ
ASGHFSVELVRGYAGFGLTLGGGRDVAGDTPLAVRGLLKDGPAQRCGRLEVGDLVLHING
ESTQGLTHAQAVERIRAGGPQLHLVIRRPLETHPGKPRGVGEPRKGVVPSWPDRSPDPGG
PEVTGSRSSSTSLVQHPPSRTTLKKTRGSPEPSPEAAADGPTVSPPERRAEDPNDQIPGS
PGPWLVPSEERLSRALGVRGAAQLAQEMAAGRRRH
>tr|A0A024QZ33|A0A024QZ33_HUMAN Coiled-coil domain containing 55, isoform CRA_a ...
MKQTKLEIQKALAEDATVYEYDSIYDEMQKKKEENNPKLLLGKDRKPKYIHNLLKAVEIR
KKEQEKRMEKKIQREREMEKGEFDDKEAFVTSAYKKKLQERAEEEEREKRAAALEACLDV
TKQKDLSGFYRHLLNQAVGEEEVPKCSFREARSGIKEEKSRGFSNEVSSKNRIPQEKCIL
QTDVKVEENPDADSDFDAKSSADDEIEETRVNCRREKVIETPENDFKHHRSQNHSRSPSE
ERGHSTRHHTKGSRTSRGHEKREDQHQQKQSRDQENHYTDRDYRKERDSHRHREASHRDS
HWKRHEQEDKPRARDQRERSDRVWKREKDREKYSQREQERDRQQNDQNRPSEKGEKEEKS
KAKEEHMKVRKERYENNDKYRDREKREVGVQSSERNQDRKESSPNSRAKDKFLDQERSNK
MRNMAKDKERNQEKPSNSESSLGAKHRLTEEGQEKGKEQERPPEAVSKFAKRNNEETVMS
ARDRYLARQMARVNAKTYIEKEDD
>tr|A0A024QZ42|A0A024QZ42_HUMAN HCG1985580, isoform CRA_c OS=Homo sapiens GN=PDCD6 PE=1 SV=1
MFDRENKAGVNFSEFTGVWKYITDWQNVFRTYDRDNSGMIDKNELKQALSGFGYRLSDQF
HDILIRKFDRQGRGQIAFDDFIQGCIVLQRLTDIFRRYDTDQDGWIQVSYEQYLSMVFSI
V
...
```

### 1B. Training peptides with PTM sites marked in lower case letters.

This provides the training information for the model, these peptides will be later mapped to the training sequences.

```
FyYEILNsPEKACSL
EAIAELDtLNEESYK
EESYKDStLIMQLLR
QLLRDNLtLWtsENQ
RDNLtLWtsENQGDE
DNLtLWtsENQGDEG
AGMDVELtVEERNLL
VEERNLLsVAyKNVI
RNLLsVAyKNVIGAR
RASWRIIssIEQKEE
ASWRIIssIEQKEEN
EYRQMVEtELKLICC
GESKVFYyKMKGDYH
MKGDYHRyLAEFATG
RKEAAENsLVAYKAA
ALNFSVFyYEILNSP
DAIAELDtLsEEsyK
IAELDtLsEEsyKDS
LDtLsEEsyKDStLI
DtLsEEsyKDStLIM
...
```

## Step 2. Three modules of PTMtopographer & Machine learning algorithm scripts

* Module for feature data generation.
* Module for prediction summarisation including calculation of false discovery rates (gDER, pDER).
* Module for annotating additional information to the protein and site-level output.

### 2A. Feature generation

This module generates feature data for candidate windows and decoy windows for the target protein sequences user provides. Suppose that we wish to train the model using protein set A(first_k) and predict PTM sites on protein set B(second_k), the user needs to run this script for both A \& B and produce the feature data for candidate windows in set A, and feature data for candidate and decoy windows in set B.

#### Input parameter file

The following items must be specified in the input parameter file for this module:


```
>length_flank=7                                                       
# number of flanking amino acids around the target site

>fasta_name=../tutorial/k_map_f_fasta.tsv           
# 1A

>pred_site=k                                                  
# target site to be predicted

>train_info=../tutorial/ubi_pep.tsv   
# 1B

>type=first_k
# or second_k
```

#### Command line

```
./program_feature_generation input_feature_generation.tsv

```

#### Output
Execution of this module on each data set (training/test) generates a number of output files for later use:

```
can_sites_position_first_k.tsv
can_sites_properties_first_k.tsv
can_sites_prots_first_k.tsv
can_sites_states_first_k.tsv
can_svm_input_first_k.tsv
decoy_sites_properties_first_k.tsv
decoy_sites_prots_first_k.tsv
decoy_svm_input_first_k.tsv
head_annotations_first_k.tsv
protein_head_annotations_first_k.tsv
my_proteins_first_k.tsv

```
Users can put these outputs to the folder  `tutorial/k_f_feature_generation_output`

Similarly, running on set B will produce the following files:

```
can_sites_position_second_k.tsv
can_sites_properties_second_k.tsv
can_sites_prots_second_k.tsv
can_sites_states_second_k.tsv
can_svm_input_second_k.tsv
decoy_sites_properties_second_k.tsv
decoy_sites_prots_second_k.tsv
decoy_svm_input_second_k.tsv
head_annotations_second_k.tsv
protein_head_annotations_second_k.tsv
my_proteins_second_k.tsv

```
Users can put these outputs to the folder  `tutorial/k_s_feature_generation_output`


### Machine learning algorithms

Before we proceed to the second module of PTMtopograher, we need to train the prediction algorithm and calculate prediction scores in the test data using the feature data generated from the first module. Here we provide the detailed steps for random forest (R script) and libsvm (C++) package.

#### Random forest

For random forest, proceed to R session.

```

###set your working directory 

setwd("PTMtopographer/")

source("randomforest.R")


training_candidate_feature="tutorial/k_f_feature_generation_output/can_sites_properties_first_k.tsv"
training_states="tutorial/k_f_feature_generation_output/can_sites_states_first_k.tsv"
test_candidate_feature="tutorial/k_s_feature_generation_output/can_sites_properties_second_k.tsv"
test_decoy_feature="tutorial/k_s_feature_generation_output/decoy_sites_properties_second_k.tsv"
outputfile_can=“tutorial/rf_prediction_second_k_can.tsv"
outputfile_decoy=“tutorial/rf_prediction_second_k_decoy.tsv"



randomforest(training_candidate_feature, training_states, test_candidate_feature, test_decoy_feature, outputfile_can, outputfile_decoy)

```

#### SVM

For training the SVM classifier, the user needs to get the best parameters in the kernel function first.
```
 grid.py can_svm_input_first_k.tsv
```

Then this should give the best gamma and cost for later training.
```
./svm-train -g 2 -c 32 -b 1 -e 0.5 -h 0 can_svm_input_first_k.tsv model_trained_on_first_k
```

Setting -b to 1 indicates that we expect the probability score as prediction output. "model_trained_on_first_k" is the name of the model we trained. Finally, we compute prediction scores for the sites on both candidate data and decoy data in protein sequences in set B(second_k).
```
./svm-predict -b 1 can_svm_input_second_k.tsv model_trained_on_first_k prediction_can_first_k_on_second_k
./svm-predict -b 1 decoy_svm_input_second_k.tsv model_trained_on_first_k prediction_decoy_first_k_on_second_k
```







### 2B. Prediction summary

The second module processes the output from the machine learning algorithm, including calculation of gDER and pDER. We illustrate this for the case of random forest prediction results here. 

#### Input parameter file

```
>mycan_prob=../tutorial/rf_prediction_second_k_can.tsv

>mydecoy_prob=../tutorial/rf_prediction_second_k_decoy.tsv

>list_prots=../tutorial/k_s_feature_generation_output/my_proteins_second_k.tsv

>mycan_prots=../tutorial/k_s_feature_generation_output/can_sites_prots_second_k.tsv

>mydecoy_prots=../tutorial/k_s_feature_generation_output/decoy_sites_prots_second_k.tsv

>can_col=0
>decoy_col=0

>decoy_tag=near_decoy
```
#### Command line
```
cd bin
```


```
./program_prediction_summary input_prediction_summary.tsv
```

#### Output

```
prot_specific_der.tsv
global_der.tsv

```
Users can put these outputs to the folder  `tutorial/k_s_prediction_summary_output`. After predicting on the first set using the second set as training, put output to the folder `tutorial/k_f_prediction_summary_output`. 





### 2C. Additional annotation of prediction output
The third module appends additional site-level and protein-level information and produces the final output of the analysis. The annotation files we provide here were generated for human proteins. Users who are performing prediction analysis on other organisms, these files must be prepared accordingly. 

#### Input parameter file
```
>prot_id=program_additional_annotation_input/protid.tsv
>dom_start=program_additional_annotation_input/env_start.tsv
>dom_end=program_additional_annotation_input/env_end.tsv
>dom_name=program_additional_annotation_input/domain_name.tsv
>dom_type=program_additional_annotation_input/domain_type.tsv

>can_sites_position=../tutorial/k_s_feature_generation_output/can_sites_position_second_k.tsv

>myprots=../tutorial/k_s_feature_generation_output/my_proteins_second_k.tsv

>prot_cytoskeleton=program_additional_annotation_input/cytoskeleton.tsv
>prot_cytosol=program_additional_annotation_input/cytosol.tsv
>prot_endoplasmic=program_additional_annotation_input/endoplasmic.tsv
>prot_endosome=program_additional_annotation_input/endosome.tsv
>prot_extracellular=program_additional_annotation_input/extracellular.tsv
>prot_golgi=program_additional_annotation_input/golgi.tsv
>prot_lysosome=program_additional_annotation_input/lysosome.tsv
>prot_membrane=program_additional_annotation_input/membrane.tsv
>prot_mitochondrion=program_additional_annotation_input/mitochondrion.tsv
>prot_nucleus=program_additional_annotation_input/nucleus.tsv
>prot_peroxisome=program_additional_annotation_input/peroxisome.tsv

>rf_score=../tutorial/rf_prediction_second_k_can.tsv

>head_annotation=../tutorial/k_s_feature_generation_output/head_annotation.tsv

>global_der=../tutorial/k_s_prediction_summary_output/global_der.tsv

>prot_specific_der=../tutorial/k_s_prediction_summary_output/prot_specific_der.tsv

>prot_head_annotation=../tutorial/k_s_feature_generation_output/protein_head_annotation.tsv

```

#### Command line
```
cd bin
```



```
./program_additional_annotation input_additional_annotation.tsv

```


#### Output

```

site_annotation.tsv
protein_annotation.tsv
```



Users can put these outputs to the folder  `tutorial/k_s_additional_annotation_output`


After generating results for the first set, users can put these outputs to the folder  `tutorial/k_f_additional_annotation_output`













































