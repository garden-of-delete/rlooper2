//
// Created by Robert Stolz on 6/27/17.
//

#ifndef RLOOPER2_STRUCTURE_H
#define RLOOPER2_STRUCTURE_H
#include <string>

class Loci{
public:
    std::string chromosome;
    std::string strand;
    long int start_pos;
    long int end_pos;
    Loci(): start_pos(0), end_pos(0) {}
    Loci(std::string, std::string, long int, long int);
    /**
     * Returns the length in base pairs of the loci
     * @returns an int representing the length of the loci
     */
    int get_length();
};

class Peak{ //probably belongs in a different .h file. Will move once I implement peak analysis. -r 6/27/17
public:
    Loci position;
    float intensity;
    Peak(): intensity(0) {}
    Peak(Loci, int);
};

class Structure{
public:
    Loci position;
    double free_energy;
    long double boltzmann_factor;
    double probability;
    double residual_twist;
    double residual_linking_difference;
    //operators
    bool operator<(const Structure &rhs) const { return free_energy < rhs.free_energy; } //overloaded < operator for sorting
    //constructors
    Structure(): free_energy(0.), boltzmann_factor(0.), probability(0.), residual_twist(0.), residual_linking_difference(0.) {}
    Structure(Loci, float, float, float);
};

#endif