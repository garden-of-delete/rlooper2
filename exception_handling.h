//
// Created by Robert Stolz on 6/28/17.
//
#ifndef RLOOPER2_EXCEPTION_HANDLING_H
#define RLOOPER2_EXCEPTION_HANDLING_H

#include <exception>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

class ModelException1: public exception{
    virtual const char* what() const throw()
    {
        return "The compute_structures function has been given a non-empty vector of structures.";
    }
};

class DefaultConstructorException : public exception{
private:
    string triggering_object;
    string message = "A prohibited default constructor has been used for object: ";

public:
    DefaultConstructorException(string object) {
        triggering_object = object;
    };

    virtual const char* what() const throw()
    {
        const string return_message = message + triggering_object;
        return return_message.c_str();
    }

};

class InvalidSequenceDataException : public exception{
private:
    string triggering_char;
    string message = "The following character in the input sequence is unrecognized: ";
public:
    InvalidSequenceDataException(char c){
        stringstream ss;
        ss << c;
        ss >> triggering_char;
    };

    virtual const char* what() const throw()
    {
        return (message + triggering_char).c_str();
    }
};

class EmptyGeneException: public exception{
    virtual const char* what() const throw()
    {
        return "There is an unrecognized or out of place charicter in the input sequence.";
    }
};

class UnexpectedEOFException: public exception{
    virtual const char* what() const throw()
    {
        return "End of file reached unexpectedly.";
    }
};

class UnexpectedClosedFileException : public exception{
private:
    string triggering_function;
    string message = "The following function has recieved a closed filestream: ";
public:
    UnexpectedClosedFileException(string c){
        stringstream ss;
        ss << c;
        ss >> triggering_function;
    };

    virtual const char* what() const throw()
    {
        return (message + triggering_function).c_str();
    }
};

class WindowerException: public exception{
    virtual const char* what() const throw()
    {
        return "A fatal error has occurred in the Windower class.";
    }
};

class SimulationException: public exception {
private:
    string description_message;
    string message = "The simulation has thrown an unhandled excpetion: ";
public:
    SimulationException(string c) {
        stringstream ss;
        ss << c;
        ss >> description_message;
    }
};
#endif //RLOOPER2_EXCEPTION_HANDLING_H
