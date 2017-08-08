//
// Created by Robert Stolz on 6/27/17.
//

#include "Rloop_equilibrium_model.h"
#include "biophysics.h"
//constructors
Rloop_equilibrium_model::Rloop_equilibrium_model() {

    N = 1500; //1500bp is the experimentally determined length of the (-) sc domain after the transcription machinery
    A = 1/10.4; // turns/bp
    C = 1.8; //tortional stiffness of ssDNA winding. (Could be 3.6 for ds or 1.8 for ss winding)
    T = 310;
    k = (2200 * 0.0019858775 * T) / N; //Hooke's law coefficient: (2200*ideal_gas_constant in kcal/mol*absolute_temp_in_kelvin)/N
    a = 5; //Neucleation Free Energy in Kcals (~3-10.2kCals) 5000
    sigma = -0.07; //measurement of energy upstream of replication domain
    alpha = N*sigma*A; //linking difference: topological parameter

    //energy difference for the described RNA/DNA dinucleotide pairs (RNA/DNA - DNA/DNA) in kcal/mol as described in Huppert 2008
    //note, energies are for the RNA strand in the 5' to 3' direction
    rGG_dCC = -0.36;
    rGC_dCG = -0.16;
    rGA_dCT = -0.1;
    rGU_dCA = -0.06;
    rCG_dGC = 0.97;
    rCC_dGG = 0.34;
    rCA_dGT = 0.45;
    rCU_dGA = 0.38;
    rAG_dTC = -0.12;
    rAC_dTG = -0.16;
    rAA_dTT = 0.6;
    rAU_dTA = -0.12;
    rUG_dAC = .45;
    rUC_dAG = .5;
    rUA_dAT = .28;
    rUU_dAA = .8;
}

//getters
int Rloop_equilibrium_model::getMinimum_loop_length() const {
    return minimum_loop_length;
}

int Rloop_equilibrium_model::getN() const {
    return N;
}

double Rloop_equilibrium_model::getA() const {
    return A;
}

double Rloop_equilibrium_model::getC() const {
    return C;
}

double Rloop_equilibrium_model::getK() const {
    return k;
}

double Rloop_equilibrium_model::geta() const {
    return a;
}

double Rloop_equilibrium_model::getSigma() const {
    return sigma;
}

double Rloop_equilibrium_model::getAlpha() const {
    return alpha;
}

void Rloop_equilibrium_model::setMinimum_loop_length(int minimum_loop_length) {
    Rloop_equilibrium_model::minimum_loop_length = minimum_loop_length;
}

void Rloop_equilibrium_model::setN(int N) {
    Rloop_equilibrium_model::N = N;
}

void Rloop_equilibrium_model::setA(double A) {
    Rloop_equilibrium_model::A = A;
}

void Rloop_equilibrium_model::setC(double C) {
    Rloop_equilibrium_model::C = C;
}

void Rloop_equilibrium_model::setK(double k) {
    Rloop_equilibrium_model::k = k;
}

void Rloop_equilibrium_model::seta(double a) {
    Rloop_equilibrium_model::a = a;
}

void Rloop_equilibrium_model::setSigma(double sigma) {
    Rloop_equilibrium_model::sigma = sigma;
}

void Rloop_equilibrium_model::setAlpha(double alpha) {
    Rloop_equilibrium_model::alpha = alpha;
}

double Rloop_equilibrium_model::getT() const {
    return T;
}

void Rloop_equilibrium_model::setT(double T) {
    Rloop_equilibrium_model::T = T;
}

double Rloop_equilibrium_model::step_forward_bps(const vector<char>::iterator& first, const vector<char>::iterator& second){
    char b_0 = *first;
    char b_1 = *second;
    double I_0;

    I_0 = compute_bps_interval(b_0, b_1); //need to verify this is the correct treatment of each interval

    return I_0;
}

double Rloop_equilibrium_model::compute_bps_interval(const char &first, const char &second){
    if (first == 'C'){ //C
        if (second == 'C') //CC
            return rGG_dCC;
        else if (second == 'G') //CG
            return rCG_dGC;
        else if (second == 'T') //CT
            return rGA_dCT;
        else //CA
            return rGU_dCA;
    }
    else if (first == 'G'){ //G
        if (second == 'C') //GC
            return rCG_dGC;
        else if (second == 'G') //GG
            return rCC_dGG;
        else if (second == 'T') //GT
            return rCA_dGT;
        else //GA
            return rCU_dGA;
    }
    else if (first == 'T'){ //T
        if (second == 'C') //TC
            return rAG_dTC;
        else if (second == 'G') //TG
            return rAC_dTG;
        else if (second == 'T') //TT
            return rAA_dTT;
        else //TA
            return rAU_dTA;
    }
    else{ //A
        if (second == 'C') //AC
            return rUG_dAC;
        else if (second == 'G') //AG
            return rUC_dAG;
        else if (second == 'T') //AT
            return rUA_dAT;
        else //AA
            return rUU_dAA;
    }
}

void Rloop_equilibrium_model::compute_structure(const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure){
    std::vector<char>::iterator b_0;
    //get boundaries of the sequence for this structure
    long int m = std::distance(start,stop);
    //compute the superhelicity term
    structure.free_energy = (2 * pow(M_PI, 2)*C*k*pow((alpha + m*A), 2)) / (4 * pow(M_PI, 2) *C + k*m);
    //compute the base-pairing energy in a loop over the sequence in the boundary
    for (b_0=start; b_0 < stop-1; b_0++){
        structure.free_energy += step_forward_bps(b_0,b_0+1);
        structure.boltzmann_factor = compute_boltzmann_factor(structure.free_energy,T);
    }
}

long double Rloop_equilibrium_model::ground_state_factor(){
    return compute_boltzmann_factor(((k*pow(alpha, 2)) / 2) - a,T);
}
