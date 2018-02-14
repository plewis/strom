bin_PROGRAMS = strom
strom_SOURCES = main.cpp strom.hpp node.hpp tree.hpp split.hpp tree_manip.hpp data.hpp gtr_model.hpp lot.hpp likelihood.hpp tree_summary.hpp xstrom.hpp
strom_CPPFLAGS = -std=c++11 -Wall \
                 -I$(HOME)/include/libhmsbeagle-1 \
                 -I$(HOME)/include \
                 -I$(HOME)/Documents/libraries/boost_1_66_0 \
                 -I$(HOME)/Documents/libraries/eigen-eigen-5a0156e40feb
strom_LDADD = -L$(HOME)/lib -lhmsbeagle \
              -L$(HOME)/lib/ncl -lncl \
              -L$(HOME)/Documents/libraries/boost_1_66_0/stage/lib -lboost_program_options
