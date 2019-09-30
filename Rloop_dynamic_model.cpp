//
// Created by Robert Stolz on 7/21/19.
//

#include "Rloop_dynamic_model.h"

//
// Created by Robert Stolz on 6/27/17.
//

//constructors
Rloop_dynamic_model::Rloop_dynamic_model() {

    N = 1500; //1500bp is the experimentally determined length of the (-) sc domain after the transcription machinery
    A = 1/10.4; // turns/bp
    C = 1.8; //tortional stiffness of ssDNA winding. (Could be 3.6 for ds or 1.8 for ss winding)
    T = 310;
    k = (2200 * 0.0019858775 * T) / N; //Hooke's law coefficient: (2200*ideal_gas_constant in kcal/mol*absolute_temp_in_kelvin)/N
    a = 10; //Neucleation Free Energy in Kcals (~3-10.2kCals) 5000
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
    homopolymer_override = false;
    unconstrained = false;
    override_energy = 0.0;
}

int Rloop_dynamic_model::getN() const {
    return N;
}

double Rloop_dynamic_model::getA() const {
    return A;
}

double Rloop_dynamic_model::getC() const {
    return C;
}

double Rloop_dynamic_model::getK() const {
    return k;
}

double Rloop_dynamic_model::geta() const {
    return a;
}

double Rloop_dynamic_model::getSigma() const {
    return sigma;
}

double Rloop_dynamic_model::getAlpha() const {
    return alpha;
}

double Rloop_dynamic_model::getAlphaTotal() const {
    return alpha_total;
}

int Rloop_dynamic_model::getCurrentPos() const {
    return current_pos;
}

void Rloop_dynamic_model::setN(int N) {
    Rloop_dynamic_model::N = N;
    setAlpha(N*sigma*A);
    k = (2200 * 0.0019858775 * T) / N;
}

void Rloop_dynamic_model::setA(double A) {
    Rloop_dynamic_model::A = A;
    setAlpha(N*sigma*A);
}

void Rloop_dynamic_model::setC(double C) {
    Rloop_dynamic_model::C = C;
}

void Rloop_dynamic_model::setK(double k) {
    Rloop_dynamic_model::k = k;
}

void Rloop_dynamic_model::seta(double a) {
    Rloop_dynamic_model::a = a;
}

void Rloop_dynamic_model::set_superhelicity(double sigma) {
    Rloop_dynamic_model::sigma = sigma;
    setAlpha(N*sigma*A);
}

void Rloop_dynamic_model::set_unconstrained(bool value){
    unconstrained = value;
}

void Rloop_dynamic_model::setAlpha(double alpha) {
    Rloop_dynamic_model::alpha = alpha;
}

void Rloop_dynamic_model::setAlphaTotal(double alpha) {
    Rloop_dynamic_model::alpha_total = alpha;
}

double Rloop_dynamic_model::getT() const {
    return T;
}

void Rloop_dynamic_model::setT(double T) {
    Rloop_dynamic_model::T = T;
    k = (2200 * 0.0019858775 * T) / N;
}

int Rloop_dynamic_model::find_distance(vector<char>& sequence, const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure){
    std::vector<char>::iterator b_0;
    int n=1;
    for (b_0=start; b_0 != stop;n++) {
        if (b_0 == sequence.end()-1) { //boundary condition
            b_0 = sequence.begin();
        }
        else{
            b_0++;
        }
    }
    return n;
}

double Rloop_dynamic_model::step_forward_bps(const vector<char>::iterator& first, const vector<char>::iterator& second){
    char b_0 = *first;
    char b_1 = *second;
    double I_0;

    I_0 = compute_bps_interval(b_0, b_1); //need to verify this is the correct treatment of each interval

    return I_0;
}

void set_bp_energy_override();

double Rloop_dynamic_model::compute_bps_interval(const char &first, const char &second){
    if (homopolymer_override == true){
        return override_energy;
    }
    if (first == 'C'){ //C
        if (second == 'C') //CC
            return rGG_dCC;
        else if (second == 'G') //CG
            return rGC_dCG;
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

void Rloop_dynamic_model::set_bp_energy_override(double energy){
    homopolymer_override = true;
    override_energy = energy;
}

void Rloop_dynamic_model::reset_model() {
    current_pos = 0;
    current_bp_energy = 0;
    current_junction_energy = 0;
    current_superhelical_energy = 0;
    proposed_bp_energy = 0;
    proposed_junction_energy = 0;
    proposed_superhelical_energy = 0;
    partition_function = 0;
    n_rloops = 0;
    n_rloop_bases = 0;
    current_pos = window_size-1;
}

void Rloop_dynamic_model::step_forward_initiation() {
    //evaluate current window
    proposed_superhelical_energy = 0;
    proposed_junction_energy = 0;
    partition_function = 0;
    for (int i=current_pos-window_size+1; i < current_pos; i++){
        int m=0;
        proposed_bp_energy = 0;
        for (int j=i+1; j <= current_pos; j++) { //sum up the boltzmann factors for the proposed energy of every rloop in the initiation window
            proposed_bp_energy += compute_bps_interval(sequence[j-1], sequence[j]);
            m = j-i+1+n_rloop_bases;
            proposed_superhelical_energy =
                    (2 * pow(M_PI, 2) * C * k * pow((alpha_total + m * A), 2)) / (4 * pow(M_PI, 2) * C + k * m);
            proposed_junction_energy = current_junction_energy + .5 * a;
            partition_function += compute_boltzmann_factor(proposed_bp_energy + proposed_superhelical_energy + proposed_junction_energy,T);
        }
    }
    partition_function += ground_state_factor();
    double urn  = ((double)rand()/(double)RAND_MAX); // random number chosen uniformly on [0,1]
    long double test = (1-ground_state_factor()/partition_function); //p(all r-loop states)
    if (urn < test){ //urn < P(all r-loop states within the initiation window
        in_rloop = true;
        n_rloops++;
        n_rloop_bases += window_size;
        //compute current energy values
        for (int i=current_pos-window_size+2; i < current_pos; i++){
            current_bp_energy += compute_bps_interval(sequence[i], sequence[i+1]);
        }
        current_junction_energy = .5 * a;
        int m = window_size;
        current_superhelical_energy =
                (2 * pow(M_PI, 2) * C * k * pow((alpha + m * A), 2)) / (4 * pow(M_PI, 2) * C + k * m);
        cout << "init " << current_pos + 2 - window_size<< endl;
        write_buffer << "init " << current_pos + 2 - window_size<< endl;
        current_loci.start_pos = current_pos + 2 - window_size;
        //print_topological_state();
    }
    else {
        current_pos += initiation_step_size; //advance the current position
    }
}

bool Rloop_dynamic_model::step_forward_elongation() {
    proposed_bp_energy = current_bp_energy;
    current_pos += elongation_step_size;
    if (current_pos >= sequence.size()-1){ //if at the end of the sequence
        return false;}
    for (int i=current_pos-elongation_step_size; i < current_pos; i++){
        proposed_bp_energy += compute_bps_interval(sequence[i],sequence[i+1]);
    }
    int m = n_rloop_bases + elongation_step_size;
    proposed_superhelical_energy =
            (2 * pow(M_PI, 2) * C * k * pow((alpha_total+ m * A), 2)) / (4 * pow(M_PI, 2) * C + k * m); // QUESTION: is this the correct term to be using here?
    current_superhelical_energy =
            (2 * pow(M_PI, 2) * C * k * pow((alpha_total + n_rloop_bases * A), 2)) / (4 * pow(M_PI, 2) * C + k * n_rloop_bases);
    proposed_junction_energy = current_junction_energy+0.5*a; //propose to add a new junction
    //proposed_junction_energy = current_junction_energy;
    long double proposed_total_energy_bf = compute_boltzmann_factor(proposed_bp_energy+proposed_superhelical_energy+current_junction_energy,T);
    partition_function += proposed_total_energy_bf;
    partition_function += compute_boltzmann_factor(current_bp_energy+current_superhelical_energy+proposed_junction_energy*a,T); //bf for the no extend state
    double urn  = ((double)rand()/(double)RAND_MAX); // random number chosen uniformly on [0,1]
    long double test = (proposed_total_energy_bf/partition_function); //p(all r-loop states)
    if (urn < test){ //if urn < P(elongation)
        n_rloop_bases += elongation_step_size;
        current_superhelical_energy = proposed_superhelical_energy;
        current_bp_energy = proposed_bp_energy;
    }
    else{ //termination
        in_rloop = false;
        current_junction_energy += .5*a;
        current_loci.end_pos = current_pos+1;
        rloop_structures.push_back(current_loci);
        cout << "term " << current_pos + 1 << ' ' << current_loci.get_length() << endl;
        write_buffer << "term " << current_pos + 1 << ' ' << current_loci.get_length() << endl;
        current_pos = current_pos + window_size; //move the current position forward to set up the initiation window
        //print_topological_state();
        return false;
    }
    return true;
}

void Rloop_dynamic_model::print_topological_state(){
    /*
    double ambient_linking_difference = N*A*sigma;
    double total_generated_linking_difference = ambient_linking_difference +
                                        (getCurrentPos()*transcriptional_superhelicity*getA());
    double residual_linking_difference = ((4*pow(pi,2)*C) / (4*pow(pi,2)*C+k*n_rloop_bases)) * (total_generated_linking_difference+n_rloop_bases*A);
    cout << "current_pos: " << current_pos+1 << endl; //convert index to bp position
    cout << "Basal linking difference: " << ambient_linking_difference << endl;
    cout << "Total generated linking difference: " << total_generated_linking_difference << endl;
    cout << "Residual superhelicity from " << n_rloops << " R-loops: " << residual_linking_difference << endl;*/

    /*
    cout << "current_pos: " << current_pos+1 << endl; //convert index to bp position
    cout << "Basal linking difference: " << ambient_linking_difference << endl;
    cout << "Total generated linking difference: " << alpha_total << endl;
    cout << "Current avaliable linking difference after "<< n_rloops << " R-loops: " << alpha << endl;
     */
    write_buffer << current_pos+1<<' '<<ambient_linking_difference<<' '<<alpha_total<<' '<<n_rloops<<' '<<alpha<< endl; //convert index to bp position
}

void Rloop_dynamic_model::compute_structure(vector<char>& sequence, const std::vector<char>::iterator &start, const std::vector<char>::iterator &stop, Structure& structure){
    std::vector<char>::iterator b_0,b_1;
    //get boundaries of the sequence for this structure
    long int m = find_distance(sequence,start,stop,structure); //need to make boundary aware, draw this value from windower
    //compute the superhelicity term
    if (!unconstrained) {
        structure.free_energy = (2 * pow(M_PI, 2) * C * k * pow((alpha + m * A), 2)) / (4 * pow(M_PI, 2) * C + k * m);
    }
    //compute the base-pairing energy in a loop over the sequence in the boundary
    for (b_0=start; b_1 != stop; b_0++){
        if (b_0 == sequence.end()){ //if you reach the end of the sequence, go back to the beginning
            b_0 = sequence.begin();
            b_1 = b_0+1;
        }
        else if (b_0 == sequence.end()-1){
            b_1 = sequence.begin();
        }
        else{
            b_1 = b_0+1;
        }
        structure.bp_energy += step_forward_bps(b_0,b_1);
    }
    structure.free_energy += structure.bp_energy;
    structure.boltzmann_factor = compute_boltzmann_factor(structure.free_energy,T);
}

void Rloop_dynamic_model::compute_external_structure(Structure& structure, Structure& rloop, Peak& external){
    std::vector<char>::iterator b_0,b_1;
    //get boundaries of the sequence for this structure
    long int m = rloop.position.get_length();
    //compute the superhelicity term
    if (!unconstrained) {
        structure.free_energy = (2 * pow(M_PI, 2) * C * k * pow((alpha + m * A), 2)) / (4 * pow(M_PI, 2) * C + k * (m-external.position.get_length())) + external.intensity;
    }
    structure.free_energy += rloop.bp_energy;
    structure.boltzmann_factor = compute_boltzmann_factor(structure.free_energy,T);

}

void Rloop_dynamic_model::compute_residuals(Structure& structure){
    structure.residual_linking_difference = ((4*pow(pi,2)*C) / (4*pow(pi,2)*C+k*structure.position.get_length()-structure.external_length)) * (alpha+structure.position.get_length()*A);
    structure.residual_twist = ((2*pi*k) / (4*pow(pi,2)*C+k*structure.position.get_length()-structure.external_length)) * (alpha+structure.position.get_length()*A);
}

double Rloop_dynamic_model::compute_residual_lk_dynamic(){
    //structure.residual_linking_difference = ((4*pow(pi,2)*C) / (4*pow(pi,2)*C+k*structure.position.get_length()-structure.external_length)) * (alpha+structure.position.get_length()*A);
    //structure.residual_twist = ((2*pi*k) / (4*pow(pi,2)*C+k*structure.position.get_length()-structure.external_length)) * (alpha+structure.position.get_length()*A);
    return ((4*pow(pi,2)*C) / (4*pow(pi,2)*C+k*n_rloop_bases)) * (alpha_total+n_rloop_bases*A);

}

void Rloop_dynamic_model::ground_state_residuals(double &twist, double &writhe){
    writhe = alpha;
    twist = k/(2*pi)*C*alpha;
}

long double Rloop_dynamic_model::ground_state_factor(){
    return compute_boltzmann_factor(((k*pow(alpha, 2)) / 2) - a,T);
}

long double Rloop_dynamic_model::ground_state_energy(){
    return ((k*pow(alpha, 2)) / 2) - a;
}


//automatically generated getters and setters (Refactor later) vvv
int Rloop_dynamic_model::getWindow_size() const {
    return window_size;
}

void Rloop_dynamic_model::setWindow_size(int window_size) {
    Rloop_dynamic_model::window_size = window_size;
}

int Rloop_dynamic_model::getN_rloops() const {
    return n_rloops;
}

void Rloop_dynamic_model::setN_rloops(int n_rloops) {
    Rloop_dynamic_model::n_rloops = n_rloops;
}

int Rloop_dynamic_model::getInitiation_step_size() const {
    return initiation_step_size;
}

void Rloop_dynamic_model::setInitiation_step_size(int initiation_step_size) {
    Rloop_dynamic_model::initiation_step_size = initiation_step_size;
}

int Rloop_dynamic_model::getElongation_step_size() const {
    return elongation_step_size;
}

void Rloop_dynamic_model::setElongation_step_size(int elongation_step_size) {
    Rloop_dynamic_model::elongation_step_size = elongation_step_size;
}

double Rloop_dynamic_model::getTranscriptional_superhelicity() const {
    return transcriptional_superhelicity;
}

void Rloop_dynamic_model::setTranscriptional_superhelicity(double transcriptional_superhelicity) {
    Rloop_dynamic_model::transcriptional_superhelicity = transcriptional_superhelicity;
}

int Rloop_dynamic_model::getN_rloop_bases() const {
    return n_rloop_bases;
}

void Rloop_dynamic_model::setN_rloop_bases(int n_rloop_bases) {
    Rloop_dynamic_model::n_rloop_bases = n_rloop_bases;
}
