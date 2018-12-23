//
// Created by Robert Stolz on 6/27/17.
//
#include "gene.h"
#include <vector>
#include "model.h"

#ifndef RLOOPER2_RLOOP_EQUILIBRIUM_MODEL_H
#define RLOOPER2_RLOOP_EQUILIBRIUM_MODEL_H

class Rloop_equilibrium_model: public Model{
protected:
    //model parameters
    //equilibrium energetics parameters
    int N; //experimentally determined length of the (-) sc domain adter the transcription machinery
    double A; //turns/bp of the B-form double helix
    double C; //tortional stiffness of ssDNA winding. (Could be 3.6 for ds or 1.8 for ss winding)
    double k; //Hooke's law coefficient: (2200*ideal_gas_constant in kcal/mol*absolute_temp_in_kelvin)/N
    double T; //temparature
    double a; //Neucleation Free Energy (junction energy) in Kcals (~3-10.2kCals) 5000
    double sigma; //measurement of energy upstream of replication domain. Moved to up the Model object.
    double alpha; //linking difference: topological parameter. refers to the percent of twists difference.
    //energy difference for the described RNA/DNA dinucleotide pairs (RNA/DNA - DNA/DNA) in kcal/mol as described in Huppert et. al. 2008
    //note, energies are for the RNA (non-template/sense DNA) strand in the 5' to 3' direction
    double rGG_dCC;
    double rGC_dCG;
    double rGA_dCT;
    double rGU_dCA;
    double rCG_dGC;
    double rCC_dGG;
    double rCA_dGT;
    double rCU_dGA;
    double rAG_dTC;
    double rAC_dTG;
    double rAA_dTT;
    double rAU_dTA;
    double rUG_dAC;
    double rUC_dAG;
    double rUA_dAT;
    double rUU_dAA;
    bool homopolymer_override;
    bool unconstrained;
    double override_energy;

public:
    //constructors
    Rloop_equilibrium_model();
    //need a special constructor that lets your specify some or all of these parameters

    int getN() const;
    double getA() const;
    double getC() const;
    double getK() const;
    double getT() const;
    double geta() const;
    double getSigma() const;
    double getAlpha() const;
    void setN(int N);
    void setA(double A);
    void setC(double C);
    void setK(double k);
    void setT(double T);
    void seta(double a);
    void set_superhelicity(double sigma);
    void set_unconstrained(bool value);
    void setAlpha(double alpha);
    void set_bp_energy_override(double energy);

    //member functions
    int find_distance(vector<char>& sequence,const vector<char>::iterator& first, const vector<char>::iterator& second, Structure& structure);
    double step_forward_bps(const vector<char>::iterator& first, const vector<char>::iterator& second);
    double compute_bps_interval(const char &first, const char &second);
    void compute_structure(vector<char>& sequence, const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure);
    void compute_residuals(Structure& structure);
    void ground_state_residuals(double& twist, double& writhe);
    long double ground_state_factor();
    long double ground_state_energy();
};
#endif //RLOOPER2_MODEL_H
