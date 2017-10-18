//
// Created by Robert Stolz on 8/8/17.
//

#ifndef RLOOPER2_LUMBERJACK_H
#define RLOOPER2_LUMBERJACK_H

#include <string>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "exception_handling.h"

using namespace std;

//logging framework
class Lumberjack {
private:
    int logging_level;
    ofstream logfile;
    string get_time();

public:
    Lumberjack(int _logging_level);
    ~Lumberjack();
    void log_error(string message);
    void log_status(string message);
    void log_debug(string message);
};


#endif //RLOOPER2_LUMBERJACK_H
