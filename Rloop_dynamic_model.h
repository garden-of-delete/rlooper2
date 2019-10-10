//
// Created by Robert Stolz on 7/21/19.
//

#ifndef RLOOPER2_RLOOP_DYNAMIC_MODEL_H
#define RLOOPER2_RLOOP_DYNAMIC_MODEL_H

#include "model.h"
#include <vector>


class Rloop_dynamic_model : public Model{
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
    double alpha_total; //total generated linking difference
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

    //dynamic model specific parameters
    double current_bp_energy;
    double current_junction_energy;
    double current_superhelical_energy;
    double proposed_bp_energy;
    double proposed_junction_energy;
    double proposed_superhelical_energy;
    double partition_function;
    int n_rloops;
    int n_rloop_bases;

protected:
    int n_simulations;
    int current_pos; //The INDEX value of the current "polymerase" position. current_pos + 1 is the base-pair number.
    int window_size;
    int initiation_step_size;
    int elongation_step_size;
    double transcriptional_superhelicity;
    Loci current_loci; //used to store the coordinates of the current rloop structure.

public:
    //constructors
    Rloop_dynamic_model();
    //need a special constructor that lets your specify some or all of these parameters

    bool in_rloop;
    vector<char> sequence;
    vector<Loci> rloop_structures;
    double ambient_linking_difference;
    stringstream write_buffer;

    //getters and setters
    int getN() const;
    double getA() const;
    double getC() const;
    double getK() const;
    double getT() const;
    double geta() const;
    double getSigma() const;
    double getAlpha() const;
    double getAlphaTotal() const;
    int getCurrentPos() const;
    int getNSimulations() const;
    void setN(int N);
    void setA(double A);
    void setC(double C);
    void setK(double k);
    void setT(double T);
    void seta(double a);
    void set_superhelicity(double sigma);
    void set_unconstrained(bool value);
    void setAlpha(double alpha);
    void setAlphaTotal(double alpha);
    void set_bp_energy_override(double energy);
    void setNSimulations(int n);

    //member functions
    void reset_model();
    void step_forward_initiation();
    bool step_forward_elongation();
    void print_topological_state();

    int find_distance(vector<char>& sequence,const vector<char>::iterator& first, const vector<char>::iterator& second, Structure& structure);
    double step_forward_bps(const vector<char>::iterator& first, const vector<char>::iterator& second);
    double compute_bps_interval(const char &first, const char &second);
    void compute_structure(vector<char>& sequence, const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure);
    void compute_external_structure(Structure& structure, Structure& rloop, Peak& external);
    void compute_residuals(Structure& structure);
    double compute_residual_lk_dynamic();
    void ground_state_residuals(double& twist, double& writhe);
    long double ground_state_factor();
    long double ground_state_energy();

    //automatically generated getters and setters vvv (refactor later)
    int getWindow_size() const;
    void setWindow_size(int window_size);
    int getN_rloops() const;
    void setN_rloops(int n_rloops);
    int getInitiation_step_size() const;
    void setInitiation_step_size(int initiation_step_size);
    int getElongation_step_size() const;
    void setElongation_step_size(int elongation_step_size);
    double getTranscriptional_superhelicity() const;
    void setTranscriptional_superhelicity(double transcriptional_superhelicity);
    int getN_rloop_bases() const;
    void setN_rloop_bases(int n_rloop_bases);
};


#endif //RLOOPER2_RLOOP_DYNAMIC_MODEL_H
