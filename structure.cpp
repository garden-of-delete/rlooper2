//
// Created by Robert Stolz on 6/27/17.
//

#include "structure.h"

Loci::Loci(std::string C, std::string c, long int s, long int e): chromosome(C), strand(c), start_pos(s), end_pos(e) {}

int Loci::get_length(){
    return end_pos - start_pos;
}

Peak::Peak(Loci l, int i): position(l), intensity(i) {}

Structure::Structure(Loci l, float f, float b, float p) {
    position = l;
    free_energy = f;
    boltzmann_factor = b;
    probability = p;
}
