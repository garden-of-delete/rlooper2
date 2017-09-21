//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_MODEL_H
#define RLOOPER2_MODEL_H

#include "structure.h"
#include "windower.h"

class Model{ //abstract class
protected:
    std::vector<char>* target_sequence; //keeps track of what gene we are in
public:
    virtual void compute_structure(const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure) {};
    virtual long double ground_state_factor() {return 0.; }; //"return 0." squashes warnings when compiled
};

#endif //RLOOPER2_MODEL_H
