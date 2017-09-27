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
        return "Fatal Exception: The compute_structures function has been given a non-empty vector of structures.";
    }
};

class DefaultConstructorException : public exception{
private:
    string triggering_object;
    string message = "Fatal Exception: You have used a prohibited default constructor for object: ";

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
    string message = "Fatal Exception: The following character in the input sequence is unrecognized: ";
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
        return "Fatal Exception: There is an unrecognized or out of place charicter in the input sequence.";
    }
};

class UnexpectedEOFException: public exception{
    virtual const char* what() const throw()
    {
        return "Fatal Exception: End of file reached unexpectedly.";
    }
};

class UnexpectedClosedFileException: public exception{
    virtual const char* what() const throw()
    {
        return "Fatal Exception: Unexpected closed file encountered by the read_gene function.";
    }
};

class WindowerException: public exception{
    virtual const char* what() const throw()
    {
        return "Fatal Exception: An error has occurred in the Windower class.";
    }
};
#endif //RLOOPER2_EXCEPTION_HANDLING_H