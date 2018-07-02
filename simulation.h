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
#include "float.h"

class Simulation{
private:
    std::vector<Model*> models;
    std::vector<Gene*> genes;
    ifstream infile;
    string outfilename;
    int minlength, power_threshold; //a minimum loop length to be applied simulation-wide.
    int top; //indicates how many of the most favorable structures to output when the --top option is used
    bool reverse_flag;
    bool complement_flag;
    bool bedfile;
    bool circular_flag;
    bool residuals;
    bool auto_domain_size;
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
    void set_circular(); //this should take a boolean
    void set_residuals(bool value);
    void set_auto_domain_size(bool value);
    void set_top(int n);
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
    void simulation_B(float superhelicity, ofstream& outfile);

    /**
     * Computes expected length for the ensemble at the given superhelicity value
     * @param superhelicity
     */
    void simulation_C(float superhelicity, ofstream& outfile);

    /**
     * WIP Dynamic simulation
     */
    void simulation_D();

    /**
    * A test environmnet for debugging purposes
    */
    void sandbox();
};

#define RLOOPER2_SIMULATION_H

#endif //RLOOPER2_SIMULATION_H
