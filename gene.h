//
// Created by Robert Stolz on 6/27/17.
//

#ifndef RLOOPER2_GENE_H
#define RLOOPER2_GENE_H
#include "structure.h"
#include "exception_handling.h"
#include "model.h"
#include <vector>
#include <iostream>
#include <fstream>

class Gene{
private:
    std::string gene_name;
    std::string header;
    Loci position;
    std::vector<char> sequence;
    std::vector<std::vector<Structure>*> structures;
    //dripseq peak related
    bool contains_peaks;
    std::vector<Peak> dripseq_peaks;

    //private member functions
    /**
     * parses FASTA headers into the Gene's Loci member.
     */
    void parse_header();

public:
    Windower windower;
    //constructors
    Gene();
    //destructors
    ~Gene();
    //Gene(std::fstream* fastafile);
    //getters and setters
    string getName();
    const string &getHeader() const;
    void setHeader(const string &header);
    const Loci &getPosition() const;
    void setPosition(const Loci &position);
    const vector<char, allocator<char>> &getSequence() const;
    //void setSequence(const vector<char, allocator<char>> &sequence);
    vector<vector<Structure>*> getStructures();
    bool isContains_peaks() const;
    void setContains_peaks(bool contains_peaks);
    const vector<Peak, allocator<Peak>> &getDripseq_peaks() const;
    void setDripseq_peaks(const vector<Peak, allocator<Peak>> &dripseq_peaks);
    //member functions
    /**
     * Reads the next FASTA record from the input file. Calls parse_header.
     * @param   fastafile   An address to the open ifstream with at least one FASTA record in it.
     * @return              returns a boolean indicating if the gene was teh last FASTA record in the file.
     */
    bool read_gene(std::ifstream& fastafile);

    /**
     * useful for debugging. Prints the header and the sequence stored in the current gene.
     */
    void print_gene();

    /**
     * Takes a model and uses it to allocate and populate a vector<Structure> of structures. Pushes the result to structures.
     * @param   model   a Model object with an implemented compute_structure method.
     */
    void compute_structures(Model& model);

    /**
     * Computes the GC skew of the provided sequence.
     * @return  returns a float representing the GC skew.
     */
    float compute_GC_skew();

    /**
     * Computes the UY skew of the provided sequence.
     * @return  returns a float representing the GC skew.
     */
    float compute_UY_skew();

    float compute_GC_content();

    /**
     * complements the sequence data (A<->T, G<->C)
     */
    void complement_sequence();

    /**
     * inverts the sequence data (GATTACA <-> ACATTAG)
     */
    void invert_sequence();

    /**
     * Returns the length of the gene
     * @return  returns an int computed from the gene's location
     */
     int get_length();

    void unload();

};


#endif //RLOOPER2_GENE_H
