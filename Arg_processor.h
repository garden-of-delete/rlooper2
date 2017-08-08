//
// Created by Robert Stolz on 8/8/17.
//

#ifndef RLOOPER2_ARG_PROCESSOR_H
#define RLOOPER2_ARG_PROCESSOR_H

#import "simulation.h"

class Arg_processor { //rudimentary functionality for now

public:
    //constructors
    Arg_processor();

    //destructors

    //member functions
    void parse_args(int argc, char* argv[], Simulation& s);

};


#endif //RLOOPER2_ARG_PROCESSOR_H
