//
// Created by Robert Stolz on 8/8/17.
//

#include "Arg_processor.h"

void Arg_processor::parse_args(int argc, char *argv[], Simulation &s){
    string infilename = argv[1];
    string outfilename = argv[2];
    s.set_infile(infilename);
    s.set_outfile(outfilename);
}