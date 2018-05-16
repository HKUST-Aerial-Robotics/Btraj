/*! \class Console
    \brief Provides high-level terminal input parsing and output visualization.

    Stand-alone, standard C++ functions encapsulated in an static class
    which helps to deal with the console input/output by displaying predefined message styles or parsing options.
    Copyright (C) 2014 Javier V. Gomez, Isaac Rivero
    www.javiervgomez.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef CONSOLE_H_
#define CONSOLE_H_

#include <string>
#include <vector>

/// \todo a member method like string but with const char* is missing.

class console {

    public:
      /**
       * Looks for the position index within the command which called the program of the parameter char[] '-whatever' and returns it,
       * or -1 if not found. Its inputs are the same as for the main function.
       *
       * NOTE: This is an internal function, not recommended to be used as a stand-alone function,
       * it is used in the parseArguments functions of this class. However, it is left public
       * because it could be useful for anybody.
       *
       * @param argc number of parameters given in the command line.
       * @param argv array of char[] with the parameters given in the command line.
       * @return The index of the char searched. -1 if not found.
       */
    static int findArguments  (int argc, const char** argv, const char* argument_name);


      /**
       * Looks for the argument given in str and stores its value in val, assuming simple input.
       *
       * Example:
       * @code int idx =  console::parseArguments (argc, argv, '-t', val); @endcode
       * and the program is run with the command:
       *  $ ./test -t 1
       * The result will be idx = 2 (third position) and val = 1.
       *
       * @param argc number of parameters given in the command line.
       * @param argv array of char[] with the parameters given in the command line.
       * @param str argument to look for.
       * @param val output variable which stores the value of the argument.
       * @return The index position in which the argument value was found.
       */
    static int parseArguments (int argc, const char** argv, const char* str, std::string &val);

     /**
      * Same as the others parseArguments. Bool as output.
      @see parseArguments()
       */
    static int parseArguments (int argc, const char** argv, const char* str, bool &val);

     /**
      * Same as the others parseArguments. int as output.
      @see parseArguments()
       */
    static int parseArguments (int argc, const char** argv, const char* str, int &val);

     /**
      * Same as the others parseArguments. float as output.
      @see parseArguments()
       */
    static int parseArguments (int argc, const char** argv, const char* str, float &val);

     /**
      * Same as the others parseArguments. double as output.
      @see parseArguments()
       */
    int parseArguments (int argc, const char** argv, const char* str, double &val);

     /**
      * Same as the others parseArguments. unsigned int as output.
      @see parseArguments()
       */
    static int parseArguments (int argc, const char** argv, const char* str, unsigned int &val);

     /**
      * Same as the others parseArguments. char as output.
      @see parseArguments()
       */
    static int parseArguments (int argc, const char** argv, const char* str, char &val);


     /**
       * Looks for the argument given in str and stores its value in val, assuming multiple input.
       *
       * Example:
       * @code
       * std::vector<std::string> vals;
       * int idx = int parseArguments (argc, argv, '-t', vals); @endcode
       * and the program is run with the command:
       *  $ ./test -input param1 param2 param3
       * The result will be idx = 2 (third position) and vals = "param1" "param2" "param3".
       *
       * @param argc number of parameters given in the command line.
       * @param argv array of char[] with the parameters given in the command line.
       * @param str argument to look for.
       * @param val output variable which stores the values of the argument, using a std::vector.<>
       * @return The index position in which the argument value was found.
       */
    static int parseArguments (int argc, const char** argv, const char* str, std::vector<std::string> & vals);

     /**
      * Same as the others parseArguments. ints as output.
      @see int parseArguments (int argc, char** argv, const char* str, std::vector<std::string> & vals);
       */
    static int parseArguments (int argc, const char** argv, const char* str, std::vector<int> & vals);

     /**
       * Predefined information message output in the console.
       *
       * @param val The string to be shown.
       * */
    static void info (const std::string &val);

     /**
       * Predefined warning message output in the console.
       *
       * @param val The string to be shown.
       * */
    static void warning (const std::string &val);

     /**
       * Predefined error message output in the console.
       *
       * @param val The string to be shown.
       * */
    static void error (const std::string &val);

     /**
       * Predefined information message returned in a std::string. Useful when
       * overloading the ostream operator (<<) or to integrate it withtin other
       * messages.
       *
       * @param val The string to be shown.
       * */
    static std::string str_info (const std::string &val);

     /**
       * Predefined warning message returned in a std::string. Useful when
       * overloading the ostream operator (<<) or to integrate it withtin other
       * messages.
       *
       * @param val The string to be shown.
       * */
    static std::string str_warning (const std::string &val);

     /**
       * Predefined error message returned in a std::string. Useful when
       * overloading the ostream operator (<<) or to integrate it withtin other
       * messages.
       *
       * @param val The string to be shown.
       * */
    static std::string str_error (const std::string &val);
};

#endif /* CONSOLE_H_ */ 
