all:bin/program_feature_generation bin/program_prediction_summary bin/program_additional_annotation


LINK.o = $(LINK.cpp)
CXXFLAGS = -Wall -O2 -std=gnu++0x -fopenmp
CXX =/usr/local/bin/g++
objs0 = program_feature_generation
objs1 = program_prediction_summary 
objs2 = program_additional_annotation
objs0 := $(addsuffix .o, $(objs0))
objs0 := $(addprefix src/, $(objs0))
objs1 := $(addsuffix .o, $(objs1))
objs1 := $(addprefix src/, $(objs1))
objs2 := $(addsuffix .o, $(objs2))
objs2 := $(addprefix src/, $(objs2))
CPPFLAGS = -I. -I$(HOME)/boost/include 


bin/program_feature_generation : $(objs0)
	$(CXX) $(objs0) -fopenmp -o $@ 

bin/program_prediction_summary: $(objs1)
	$(CXX) $(objs1) -fopenmp -o $@ 

bin/program_additional_annotation: $(objs2)
	$(CXX) $(objs2) -fopenmp -o $@ 



clean:
	$(RM) bin/program_feature_generation bin/program_prediction_summary bin/program_additional_annotation $(objs) 
	$(RM) src/program_feature_generation.o bin/program_prediction_summary.o bin/program_additional_annotation.o $(objs) 
	

