emd_flow: main.cc emd_flow.o emd_flow.h
	g++ -Wall -Wextra -O2 -o emd_flow main.cc emd_flow.o -lboost_program_options

emd_flow.o: emd_flow.cc emd_flow.h
	g++ -Wall -Wextra -fPIC -O2 -c -o emd_flow.o emd_flow.cc -lemon -I lemon/include -L lemon/lib

mexfile: emd_flow.o emd_flow.h mex_wrapper.cc
	mex -output emd_flow mex_wrapper.cc emd_flow.o

all: emd_flow mexfile
