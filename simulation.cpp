//
// Created by Robert Stolz on 6/28/17.
//
#include "simulation.h"

void initialize_models(char* model){

}

Simulation::Simulation(){
    //default constructor
    minlength = 2; //default minlength
    reverse_flag = false;
    complement_flag = false;
    power_threshold = 1;
    circular_flag = false;
}

Simulation::~Simulation(){
    for(std::vector<Gene*>::iterator it = genes.begin(); std::distance(genes.begin(),it) < genes.size(); ++it){
        delete *it; //need to test this destructor
    }
}

void Simulation::set_infile(string infilename){
    infile.open(infilename, ios::in);
}

void Simulation::set_outfile(string Outfilename){
    outfilename = Outfilename;
}

void Simulation::set_minlength(int Minlength){
    minlength = Minlength;
}

void Simulation::set_bedfile(bool value){
    bedfile = value;
}

void Simulation::set_power_threshold(int Power_threshold){
    power_threshold = Power_threshold;
}

void Simulation::set_circular(){
    circular_flag = true;
}

void Simulation::reverse_input(){
    reverse_flag = true;
}

void Simulation::complement_input(){
    complement_flag = true;
}

std::vector<Model*> Simulation::get_models(){
    return models;
}

void Simulation::add_model(Model& model){
    models.push_back(&model);
}

void Simulation::compute_signal_bpprobs(Gene &gene, vector<double> *&signal){
    signal = new vector<double>(gene.get_length(), 0.0);
    //compute the r-loop involvement probability for each base
    //for each structure in the gene
    for (std::vector<Structure>::iterator it = gene.getStructures()[0]->begin();
         it < gene.getStructures()[0]->end(); ++it) {
        //for each base in the structure
        for (long int i = it->position.start_pos - gene.getPosition().start_pos;
             i < it->position.end_pos - gene.getPosition().start_pos; i++) {
            (*signal)[i] += it->probability;
        }
    }
    //if strand is -, reverse bp_probabilities
    if (gene.getPosition().strand == "-") {
        std::reverse(signal->begin(), signal->end());
    }
}

void Simulation::compute_signal_average_G(Gene &gene, vector<double> *&signal){
    signal = new vector<double>(gene.get_length(), 0.0);
    //compute the special partition function for each base-pair
    vector<double> bp_partition_functions(gene.get_length(), 0.0);
    //for each structure in the gene
    for (std::vector<Structure>::iterator it = gene.getStructures()[0]->begin();
         it < gene.getStructures()[0]->end(); ++it) {
        //for each base in the structure
        for (long int i = it->position.start_pos - gene.getPosition().start_pos;
             i < it->position.end_pos - gene.getPosition().start_pos; i++) {
            bp_partition_functions[i] += it->boltzmann_factor;
        }
    }
    //compute the r-loop involvement probability for each base (will probably be moved out of this func later)
    //for each structure in the gene
    for (std::vector<Structure>::iterator it = gene.getStructures()[0]->begin();
         it < gene.getStructures()[0]->end(); ++it) {
        //for each base in the structure
        for (long int i = it->position.start_pos - gene.getPosition().start_pos;
             i < it->position.end_pos - gene.getPosition().start_pos; i++) {
            (*signal)[i] += (it->boltzmann_factor/bp_partition_functions[i])*it->free_energy;
        }
    }
    //if strand is -, reverse signal
    if (gene.getPosition().strand == "-") {
        std::reverse(signal->begin(), signal->end());
    }
}

void Simulation::compute_signal_mfe(Gene &gene, vector<double> *&signal){
    signal = new vector<double>(gene.get_length(), 0.0);
    double current_min = FLT_MAX;
    Structure mfe;
    //for each structure in the gene
    for (std::vector<Structure>::iterator it = gene.getStructures()[0]->begin();
         it < gene.getStructures()[0]->end(); ++it) {
        if (it->free_energy < current_min){
            current_min = it->free_energy;
            mfe = *it;
        }
    }
    //record the position of the mfe into the signal
    for (long int i = mfe.position.start_pos - gene.getPosition().start_pos;
         i < mfe.position.end_pos - gene.getPosition().start_pos; i++) {
        (*signal)[i] = 1.0;
    }
    //if strand is -, reverse signal
    if (gene.getPosition().strand == "-") {
        std::reverse(signal->begin(), signal->end());
    }
}


void Simulation::call_peaks_threshold(Gene& gene, vector<double>& signal, vector<Loci>& peaks){
    //int power_threshold = 12; //needs to be made a class variable
    double minimum = 1;
    bool in_peak = false;
    long peak_start=0, peak_end=0;
    double magnitude = 0;
    Structure* temp;
    for (int i=0; i < signal.size(); i++){
        //determine lowest value in the signal
        if (signal[i] < minimum){
            minimum = signal[i];
        }
    }
    for (int i=0; i < signal.size(); i++){
        if (signal[i] > minimum*pow(10,power_threshold)){ //the signal is significant
            if (!in_peak){
                in_peak = true;
                peak_start = gene.getPosition().start_pos + i;
            }
        }
        else{ //the signal is not significant
            if (in_peak){
                in_peak = false;
                peak_end = gene.getPosition().start_pos + i;
                peaks.emplace_back(Loci(gene.getPosition().chromosome,gene.getPosition().strand, peak_start, peak_end)); //chromosome, strand, start_pos, end_pos
            }
        }
    }
}

void Simulation::cluster_k_intervals(vector<Loci> &peaks, vector<Loci> &clustered_peaks){
    long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::cout << "rng seed: " << seed << endl;
    vector<double> costs;
    vector<int> chosen_peaks;
    vector<int> clustering_tally;
    for (int i=0; i<peaks.size(); i++){
        clustering_tally.push_back(0);
    }
    int k;
    k = 5;
    for (int i=0; i < 1000; i++){
        lloyds_algorithm(peaks,chosen_peaks,k,seed);
        for (int j=0; j < chosen_peaks.size(); j++){
            clustering_tally[chosen_peaks[j]]++;
        }
        chosen_peaks.empty();
    }
    //push the most common cluster representatives onto clustered peaks
}

double Simulation::lloyds_algorithm(vector<Loci> &peaks, vector<int> &clustering, int k, unsigned seed){
    bool swaps = true;
    vector<int> medoid_indeces; //maps medoid index to actual element in the matrix
    vector<int> medoid_assignments; //the INDEX of the medoid each peak is assigned to.
    vector<vector<double>> pairwise_distance_matrix;
    double configuration_cost = 0;
    for (int i=0; i < peaks.size(); i++){ //initialize the pairwise distance matrix
        vector<double> temp;
        for (int j=0; j < peaks.size(); j++){
            temp.push_back(0);
        }
        pairwise_distance_matrix.push_back(temp);
    }
    //choose k different intervals at random as the initial medoids
        //generate k random indeces
    vector<int> shuffled;
    for (int i=0;i<peaks.size();i++){ //unshuffled medoid indeces
        shuffled.push_back(i);
        medoid_assignments.push_back(0); //all peaks are temporarily assigned to the first medoid
    }
    std::shuffle(shuffled.begin(),shuffled.end(),std::default_random_engine(seed)); //not tested, need to connect the seed
    for (int i=0;i<k;i++){
        medoid_indeces.push_back(shuffled[i]); //save the k randomly selected medoid indeces to a list
    }
    //compute the pairwise distance matrix
    for (int i=0;i < peaks.size();i++) { //for each peak
        for (int j=0; j < peaks.size(); j++) { //for each peak
            pairwise_distance_matrix[i][j] = interval_distance(peaks[i], peaks[j]);
        }
    }
    double current_cost = 0; //cost of the current clustering configuration
    //assign each interval to its closest medoid
    for (int i=0; i < peaks.size();i++){ //for each peak
        for (int j=1; j<k; j++){ //for each medoid index
            if(pairwise_distance_matrix[i][medoid_indeces[j]] < pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]]){
                medoid_assignments[i] = j;
            }
        }
    }
    //compute full configuration cost
    for (int i=0; i<medoid_assignments.size(); i++){
        configuration_cost += pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]];
    }

    while (swaps) { //Veroni descent
        swaps = false;
        //assign each interval number to its closest medoid (already done for the first iteration)
        for (int i=0; i < peaks.size();i++){
            for (int j=1; j<k; j++){
                if(pairwise_distance_matrix[i][medoid_indeces[j]] < pairwise_distance_matrix[i][medoid_assignments[i]]){
                    //configuration_cost -= pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]]; //update the configuration cost
                    medoid_assignments[i] = j; //update the medoid assignment with the index of the new medoid
                    //configuration_cost += pairwise_distance_matrix[i][medoid_indeces[medoid_assignments[i]]];
                }
            }
        }
        //for each cluster
        for (int p=0; p < k; p++){
            //test each object within the cluster as the new medoid of the cluster
            for (int i=0; i < peaks.size(); i++){
                if (medoid_assignments[i] == p && i != medoid_indeces[p]){ //if the medoid is in the currently considered group, but is not the current medoid
                    //determine swap cost
                    double costA = 0, costB=0;
                    for (int j=0; j < peaks.size(); j++){
                        if (medoid_assignments[i] == p) { //if element is in the currently considered cluster
                            costA += pairwise_distance_matrix[medoid_indeces[medoid_assignments[i]]][j]; //current configuration
                            costB += pairwise_distance_matrix[i][j]; //currently considered swap
                        }
                        if (costB < costA){ //swap would reduce the configuration cost
                            //update the configuration cost
                            configuration_cost -= costA;
                            configuration_cost += costB;
                            //update medoid_indeces
                            medoid_indeces[p] = i;
                        }
                    }
                }
            }
        }
    }
    //tally the final clustering
    clustering = medoid_indeces;
    return configuration_cost;
}

double Simulation::compute_configuration_cost(vector<vector<double>> &pairwise_distance_matrix,
                                              vector<int> medoid_indeces) {
    double configuration_cost = 0;
    for (int i=0; i<pairwise_distance_matrix.size(); i++){
        for (int j=0; j<medoid_indeces.size();j++) {
            configuration_cost += pairwise_distance_matrix[i][medoid_indeces[j]];
        }
    }
    return configuration_cost;
}

double Simulation::interval_distance(const Loci &A, const Loci &B){
    double term1 = pow((A.start_pos+A.end_pos)/2.-(B.start_pos+B.end_pos)/2.,2);
    double term2 = pow((A.end_pos-A.start_pos)/2.-(B.end_pos-B.start_pos)/2.,2)/3.;
    return term1+term2;
}

void Simulation::write_wigfile_header(ofstream& outfile, string trackname){
    //open stringstream
    std::stringstream ss;
    //compose .wig header
    //adjust browser position
    ss << "track type=wiggle_0 name=\"" << trackname << "\" visibility=full autoscale=off color=50,150,255 priority=10"
       << endl;
    outfile << ss.rdbuf();
}

void Simulation::write_wigfile(ofstream& outfile, Gene* gene, std::vector<double>* signal){
    //open stringstream
    std::stringstream ss;
    string wigfile_name = gene->getHeader().c_str();
    //compose .wig header
    string name = gene->getName();
    //adjust browser position
    ss << "browser position " << gene->getPosition().chromosome << ':' << gene->getPosition().start_pos << '-' <<
       gene->getPosition().end_pos << endl;
    ss << "fixedStep chrom=" << gene->getPosition().chromosome << " start=" << gene->getPosition().start_pos << " step=1"
       << endl;
    for (int i = 0; i < signal->size(); i++) {
        ss << (*signal)[i] << endl;
    }
    //write stringstream to file
    outfile << ss.rdbuf();
}

void Simulation::read_bedfile(ifstream &bedinput, vector<Loci> &peaks){
    Loci temp;
    long int pos;
    char buffer[1000];
    string strbuff;
    if (!bedinput.is_open()){
        //throw exception
    }
    while(bedinput.getline(buffer,1000)){
        strbuff = std::string(buffer);
        //need to deal with lines that do not contain a bed entry here

        //parse out chromosome name
        pos = strbuff.find('\t');
        temp.chromosome = strbuff.substr(0,pos); //need to handle non-numeric chromosome names as well
        strbuff = strbuff.substr(pos+1,strbuff.length());
        //parse out start position of the entry
        pos = strbuff.find('\t');
        temp.start_pos = stol(strbuff.substr(0,pos));
        strbuff = strbuff.substr(pos+1,strbuff.length());
        //parse out end position of the entry
        pos = strbuff.find('\t');
        temp.end_pos = stol(strbuff.substr(0,pos));
        strbuff = strbuff.substr(pos+1,strbuff.length());
        //discard the next two columns (may need to be made more flexible in the future)
        pos = strbuff.find('\t');
        strbuff = strbuff.substr(pos+1,strbuff.length());
        pos = strbuff.find('\t');
        strbuff = strbuff.substr(pos+1,strbuff.length());
        //parse out the strand
        pos = strbuff.find('\t');
        temp.strand = strbuff.substr(0,pos);
        strbuff = strbuff.substr(pos+1,strbuff.length());
        //save to the peaks vector
        peaks.push_back(temp);
    }
}

void Simulation::write_bedfile_header(ofstream& outfile, string trackname){
    //write bedfile
    stringstream ss;
    ss << "track name=rLooper description=\""<<trackname<<"\" useScore=1" << endl;
    outfile << ss.rdbuf();
}

void Simulation::write_bedfile(ofstream& outfile, Gene* gene, vector<Loci>& peaks){
    //write bedfile
    stringstream ss;
    string strand_name;
    int start_pos=0, end_pos=0;
    if (gene->getPosition().strand == "+"){
        strand_name = "POS";
    }
    else {
        strand_name = "NEG";
    }
    ss << "browser position " << gene->getPosition().chromosome << ':' << gene->getPosition().start_pos << '-' <<
       gene->getPosition().end_pos << endl;
    //print BED header here
    //print the peaks in BED format
    for (int i=0; i < peaks.size(); i++){
        ss << peaks[i].chromosome << '\t' << (peaks)[i].start_pos << '\t' << peaks[i].end_pos
           << '\t' << strand_name << i << '\t' << '0' << '\t' << peaks[i].strand << endl;
    }
    //write stringstream to file
    outfile << ss.rdbuf();
}

void Simulation::simulation_A(){ //some of this code might be migrated into new objects and functions in the future
    //initialize variables
    if (!infile.is_open()){
        throw UnexpectedClosedFileException("Simulation::simulation_A");
    }
    ofstream outfile1(outfilename+"_bpprob.wig",ios::out);
    ofstream outfile2(outfilename+"_avgG.wig",ios::out);
    ofstream outfile3(outfilename+"_mfe.wig",ios::out);
    ofstream outfile4(outfilename+"_bpprob.bed",ios::out);
    ofstream outfile6(outfilename+"_mfe.bed",ios::out);

    //write headers
    write_wigfile_header(outfile1,"signal1_"+outfilename);
    write_wigfile_header(outfile2,"signal2_"+outfilename);
    write_wigfile_header(outfile3,"signal3_"+outfilename);
    write_bedfile_header(outfile4,"signal1_peaks_"+outfilename);
    write_bedfile_header(outfile6,"signal3_peaks_"+outfilename);

    bool eof = false;
    if (models.size() < 1){
        //throw exception
    }
    //do while !eof
    while(eof == false) {
        //allocate new gene
        Gene *this_gene = new Gene();
        this_gene->windower.set_min_window_size(minlength);
        //read gene
        eof = this_gene->read_gene(infile);
        cout << "processing gene: " << this_gene->getName() << "...";
        //compute structures using models
        if (this_gene->getPosition().strand == "+") {
            //this_gene->complement_sequence();
        }
        else if(this_gene->getPosition().strand == "-") {
            this_gene->invert_sequence();
        }
        if (complement_flag) {
            this_gene->complement_sequence();
        }
        if (reverse_flag) {
            this_gene->invert_sequence();
        }
        if (circular_flag) {
            this_gene->compute_structures_circular(*models[0]);
        }
        else{
            this_gene->compute_structures(*models[0]);
        }

        //ensemble analysis, free energies and boltzmann factors have already been computed in compute_structures
        //compute partition function
        long double partition_function = 0;
        long double sanity_check = 0;
        for (vector<Structure>::iterator it = this_gene->getStructures()[0]->begin();
             it < this_gene->getStructures()[0]->end(); ++it){
               partition_function += it->boltzmann_factor;
        }
        partition_function += models[0]->ground_state_factor();
        //compute boltzmann weights and store in the structures
        for (vector<Structure>::iterator it = this_gene->getStructures()[0]->begin();
             it < this_gene->getStructures()[0]->end(); ++it){
            it->probability = it->boltzmann_factor/partition_function;
            sanity_check += it->boltzmann_factor/partition_function;
        }
        sanity_check += models[0]->ground_state_factor()/partition_function;
        double test = fabs(1-sanity_check);
        if (fabs(1-sanity_check) > .00001){
            throw SimulationException("Ensemble probability sum != 1"); //this throw is uncaught
        }
        //compute signals and output .wig tracks
        vector<double>* signal = NULL, *signal2 = NULL, *signal3 = NULL;
        vector<Loci> peaks;
        compute_signal_bpprobs(*this_gene,signal);
        compute_signal_average_G(*this_gene,signal2);
        compute_signal_mfe(*this_gene,signal3);
        write_wigfile(outfile1,this_gene,signal);
        write_wigfile(outfile2,this_gene,signal2);
        write_wigfile(outfile3,this_gene,signal3);
        //call peaks and write results to .bed files
        if (bedfile){
            call_peaks_threshold(*this_gene,*signal,peaks); //possible null pointer exception generated here
            //write to bedfile
            write_bedfile(outfile4,this_gene,peaks);
            peaks.clear();
            call_peaks_threshold(*this_gene,*signal3,peaks); //possible null pointer exception generated here
            //write to bedfile
            write_bedfile(outfile6,this_gene,peaks);
        }
        cout << "complete!" << endl;
        delete signal;
        //clear_sequence the sequence data from the gene to save memory
        this_gene->clear_sequence();
        this_gene->clear_structures();
        //store the gene in the genes vector
        genes.push_back(this_gene);
    }
    outfile1.close();
    outfile2.close();
    outfile3.close();
    outfile4.close();
}

//computes P(R-Loop) for the provided supercoiling value
void Simulation::simulation_B(float superhelicity, ofstream& outfile){
    if (!infile.is_open()){
        throw UnexpectedClosedFileException("Simulation::simulation_B");
    }
    if (models.size() < 1){
        //throw exception
    }
    float p_rloop = 0;
    Gene* this_gene;
    if (!genes.size()){
        this_gene = new Gene();
        this_gene->read_gene(infile);
        this_gene->windower.set_min_window_size(minlength);
        this_gene->complement_sequence();
        //this_gene->invert_sequence();
        genes.push_back(this_gene);
    }
    else{
        this_gene = genes[0];
    }
    models[0]->set_superhelicity(superhelicity); //set the superhelicity in the model to the provided value
    this_gene->clear_structures(); //saves memory
    this_gene->compute_structures(*(models[0]));
    //determine P(ground state)
    long double partition_function = 0;
    long double ground_state_factor = 0;
    int index  = this_gene->getStructures().size();
    int count = 0;
    for (vector<Structure>::iterator it = this_gene->getStructures()[index-1]->begin();
         it < this_gene->getStructures()[index-1]->end(); ++it){
        partition_function += it->boltzmann_factor;
        count++;
    }
    ground_state_factor = models[0]->ground_state_factor();
    partition_function += ground_state_factor;
    //determine P(R-Loop) as 1-P(ground state)
    p_rloop = 1 - (ground_state_factor/partition_function);
    //display result
    outfile << superhelicity << ' ' << p_rloop << endl;
}

void Simulation::simulation_C(float superhelicity, ofstream& outfile){
    if (!infile.is_open()){
        throw UnexpectedClosedFileException("Simulation::simulation_C");
    }
    if (models.size() < 1){
        //throw exception
    }
    Gene* this_gene;
    if (!genes.size()){
        this_gene = new Gene();
        this_gene->read_gene(infile);
        this_gene->windower.set_min_window_size(minlength);
        this_gene->complement_sequence();
        //this_gene->invert_sequence();
        genes.push_back(this_gene);
    }
    else{
        this_gene = genes[0];
    }
    models[0]->set_superhelicity(superhelicity); //set the superhelicity in the model to the provided value
    this_gene->compute_structures(*models[0]);
    //determine P(ground state)
    long double partition_function = 0;
    long double ground_state_factor = 0;
    unsigned long index = this_gene->getStructures().size();
    int count = 0;
    for (vector<Structure>::iterator it = this_gene->getStructures()[index-1]->begin();
         it < this_gene->getStructures()[index-1]->end(); ++it){
        partition_function += it->boltzmann_factor;
        count++;
    }
    ground_state_factor = models[0]->ground_state_factor();
    partition_function += ground_state_factor;
    //determine expected length at the given superhelicity value
    double expected_length = 0;
    for (vector<Structure>::iterator it = this_gene->getStructures()[index-1]->begin();
         it < this_gene->getStructures()[index-1]->end(); ++it){
        expected_length += (it->boltzmann_factor/partition_function)*it->position.get_length();
    }
    //display result
    outfile << superhelicity << ' ' << expected_length << endl;
}

void Simulation::sandbox() { //test/debug environment
/*
    srand(454); //needs to be an argument
    ifstream test("test.bed",ios::in);
    vector<Loci> testvector, clustered_peaks;
    read_bedfile(test,testvector);
    cluster_k_intervals(testvector,clustered_peaks);
*/
    /*
     * Craig's graph function is here. needs to be migrated elsewhere so sandbox can continue being used as a test function.
     *
     * */

    ofstream outfile(outfilename, ios::out);
    if (!infile.is_open()) {
        throw UnexpectedClosedFileException("Simulation::sandbox");
    }
    Gene geneA;
    geneA.read_gene(infile); //sequence being read in is not used for anything
    //craig's simulation
    double lower_bound = -0.5;
    double upper_bound = 1.0;
    double supercoiling = 0.0;
    double last_supercoiling = 0.0;
    vector<double> x;
    vector<double> y;
    double step_size = .001;
    double tolerance = 0.05;
    char last_direction = 'n';

    Rloop_equilibrium_model modelA;
    //modelA.setMinimum_loop_length(minlength); //not functional, needs to be removed
    geneA.windower.set_min_window_size(minlength);
    //for each base pairing energy
    for (double bp_energy = lower_bound; bp_energy <= upper_bound; bp_energy += 0.01) {
        //for each level of supercoiling
        supercoiling = last_supercoiling;
        while (true) {
            cout << "For bp_energy: " << bp_energy << ", and supercoiling: " << supercoiling << endl;
            //set supercoiling
            modelA.set_superhelicity(supercoiling);
            modelA.set_bp_energy_override(bp_energy);
            //dump previous ensemble of structures
            geneA.clear_structures();
            //compute ensemble
            geneA.compute_structures(modelA);
            //sum the probabilities of all R-loop structures
            double partition_function = 0;
            for (vector<Structure>::iterator it = geneA.getStructures()[0]->begin();
                 it != geneA.getStructures()[0]->end(); ++it) {
                partition_function += it.base()->boltzmann_factor;
            }
            //compare to the partition function
            // if the P(R-loop) is > 50% + tol
            if (partition_function / (partition_function + modelA.ground_state_factor()) < .5) {
                //adjust supercoiling down
                if (last_direction == 'u') {
                    //overshot the target and outside the error tolerance
                    step_size /= 0.5; //cut step size in half
                }
                supercoiling -= step_size;
                last_direction = 'd';
            }
                //if the {(R-loop) is < 50% - tol
                //adjust supercoiling up
                /*
                    if (last_direction == 'd'){
                        //overshot the target and outside the error tolerance
                        step_size /= 0.5; //cut step size in half
                    }
                    supercoiling += step_size;
                    last_direction = 'u';
                } *///end comment block here
            else { //hit target within tolerance
                cout << "complete" << endl;
                x.push_back(bp_energy);
                y.push_back(supercoiling);
                last_supercoiling = supercoiling; //improves search efficiency at the next bp_energy
                break;
            }
        }
    }
    //write result
    stringstream ss;
    for (int i = 0; i < x.size(); i++) {
        ss << x[i] << ' ' << y[i] << endl;
        outfile << ss.rdbuf();
    }
    outfile.close();
}