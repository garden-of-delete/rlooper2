//
// Created by Robert Stolz on 6/28/17.
//

#ifndef RLOOPER2_SIMULATION_H
#include "Rloop_equilibrium_model.h"
#include <algorithm>
#include <array>
#include <random>
#include <chrono>
#include "exception_handling.h"
#include <sstream>
#include <cstdlib>
#include <cmath>

class Simulation{
private:
    std::vector<Model*> models;
    std::vector<Gene*> genes;
    ifstream infile;
    string outfilename;
    int minlength, power_threshold; //a minimum loop length to be applied simulation-wide.
    bool reverse_flag;
    bool complement_flag;
    bool bedfile;
    //member functions
    void compute_signal_bpprobs(Gene &gene, vector<double> *&signal);
    void compute_signal_average_G(Gene &gene, vector<double> *&signal);
    void compute_signal_mfe(Gene &gene, vector<double> *&signal);
    void call_peaks_threshold(Gene& gene, vector<double>& signal, vector<Loci>& peaks);
    /**
     * writes a given signal for a given gene to the outfile in .wig format
     * @param   gene    the gene from which to derive the header in the wigfile
     * @param   signal  the signal to be written to wigfile
     *
     */
    //clustering functions
    void cluster_k_intervals(vector<Loci>& peaks, vector<Loci>& clustered_peaks);
    double lloyds_algorithm(vector<Loci>& peaks, vector<int>& clustering, int k, unsigned seed);
    /**
     * computes the distance metric between two interval numbers given in Guo et. al. 2014
     * @param A the firt interval number represented as a Loci object
     * @param B the second interval number represented as a Loci object
     * @return the distance between A and B
     */
    double compute_configuration_cost(vector<vector<double>>& pairwise_distance_matrix, vector<int> medoid_indeces);
    double interval_distance(const Loci &A, const Loci &B);
    void write_wigfile_header(ofstream& outfile, string trackname);
    void write_wigfile(ofstream& outfile, Gene* gene, std::vector<double>* signal);
    void read_bedfile(ifstream& bedinput, vector<Loci>& peaks);
    void write_bedfile_header(ofstream& outfile, string trackname);
    void write_bedfile(ofstream& outfile, Gene* gene, vector<Loci>& peaks);

public:
    Simulation();
    Simulation(int argc, char* argv[]);
    ~Simulation();
    //getters and setters
    void set_infile(string infilename);
    void set_outfile(string outfilename);
    void set_bedfile(bool value);
    void set_minlength(int Minlength);
    void set_power_threshold(int Power_threshold);
    void reverse_input();
    void complement_input();
    std::vector<Model*> get_models();
    void add_model(Model& model);

    //simulation protocols
    /**
     * Ensemble analysis on a set of genes. Supports one model.
     */
    void simulation_A();

    /**
     * Computes P(R-Loop is on the sequence) for a given superhelicity level. Supports any number of models.
     */
    void simulation_B(float superhelicity);

    /**
     * Computes expected length for the ensemble at the given superhelicity value
     * @param superhelicity
     */
    void simulation_C(float superhelicity);

    /**
    * A test environmnet for debugging purposes
    */
    void sandbox();
};

#define RLOOPER2_SIMULATION_H

#endif //RLOOPER2_SIMULATION_H
