/*
* console.cpp
* Helper class to console input-output.
* Copyright (C) 2014 Javier V. Gomez, Isaac Rivero
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include "console.h"

using namespace std;

//IMPORTANT NOTE: findArguments functions have not been tested (just one of them).
//It is possible to have bugs here with the index numbering.

// This function looks for the position index of the parameter string "-x" and return it (-1 if not found).
int console::findArguments 
(int argc, const char** argv, const char* argument_name) {
    for (int i = 1; i < argc; ++i)
        // Search for the string
        if (strcmp (argv[i], argument_name) == 0)
            return (i);

    return (-1);
}


// All these functions works in the same way: they look for the position of the string "-x" so index is the
// index of argv in which the value for "-x" is set. For example: ./test -t 1 
// index will be 2 (referring to the value 1).


////////////////////////////////////////////////////////////////////////////////
int console::parseArguments 
(int argc, const char** argv, const char* str, string &val) {
    int index = findArguments (argc, argv, str) + 1;
    if (index > 0 && index < argc)
        val = argv[index];
    return index;
}


////////////////////////////////////////////////////////////////////////////////
int console::parseArguments
(int argc, const char** argv, const char* str, bool &val) {
    int index = findArguments (argc, argv, str) + 1;
    if (index > 0 && index < argc )
        val = atoi (argv[index]) == 1;
    return index;
}
////////////////////////////////////////////////////////////////////////////////
int console::parseArguments 
(int argc, const char** argv, const char* str, double &val) {
    int index = findArguments (argc, argv, str) + 1;
    if (index > 0 && index < argc )
        val = atof (argv[index]);
    return index;
}
////////////////////////////////////////////////////////////////////////////////
int console::parseArguments 
(int argc, const char** argv, const char* str, float &val) {
    int index = findArguments (argc, argv, str) + 1;
    if (index > 0 && index < argc )
        val = static_cast<float> (atof (argv[index]));
    return index;
}
////////////////////////////////////////////////////////////////////////////////
int console::parseArguments 
(int argc, const char** argv, const char* str, int &val) {
    int index = findArguments (argc, argv, str) + 1;
    if (index > 0 && index < argc )
        val = atoi (argv[index]);
    return index;
}
////////////////////////////////////////////////////////////////////////////////
int console::parseArguments 
(int argc, const char** argv, const char* str, unsigned int &val) {
    int index = findArguments (argc, argv, str) + 1;
    if (index > 0 && index < argc )
        val = atoi (argv[index]);
    return index;
}
////////////////////////////////////////////////////////////////////////////////
int console::parseArguments 
(int argc, const char** argv, const char* str, char &val) {
    int index = findArguments (argc, argv, str) + 1;
    if (index > 0 && index < argc )
        val = argv[index][0];
    return index;
}


////////////////////////////////////////////////////////////////////////////////

int console::parseArguments 
(int argc, const char** argv, const char* str, vector<string> & vals) {
int index = findArguments (argc, argv, str);
    int i = index + 1;
    string s;
    do {
        const char *aux = argv[i];
        if (*aux == '-')
            i = argc;
        else {
            s = argv[i];
            vals.push_back(s);
            ++i;
        }
    } while (i < argc);
    return index;
}

int console::parseArguments 
(int argc, const char** argv, const char* str, vector<int> & vals) {
    int index = findArguments (argc, argv, str);
    int i = index + 1;
    int val;
    do {
        const char *aux = argv[i];
        if (*aux == '-')
            i = argc;
        else {
            val = atoi(argv[i]);
            vals.push_back(val);
            ++i;
        }
    } while (i < argc);
    return index;
}


void console::info
(const std::string &val) {
    cout << "\33[1;34m" << "[INFO] " << val << "\33[0m" << endl;
}


void console::warning
(const std::string &val) {
    cout << "\33[1;33m" << "[WARNING] " << val << "\33[0m" << endl;
}

void console::error
(const std::string &val) {
    cout << "\33[1;31m" << "[ERROR] " << val << "\33[0m" << endl;
}


string console::str_info 
(const string &val){
    stringstream ss;
    ss << "\33[1;34m" << "[INFO] " << val << "\33[0m" << endl;
    return ss.str();
}

string console::str_warning
(const std::string &val) {
    stringstream ss;
    ss << "\33[1;33m" << "[WARNING] " << val << "\33[0m" << endl;
    return ss.str();
}

string console::str_error
(const std::string &val) {
    stringstream ss;
    ss << "\33[1;31m" << "[ERROR] " << val << "\33[0m" << endl;
    return ss.str();
}
