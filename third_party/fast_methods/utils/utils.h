/*! \class utils
    \brief Provides helper code not related with specific classes.
    Copyright (C) 2014 Javier V. Gomez
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

#ifndef UTILS_H_
#define UTILS_H_

#include <limits>

class utils {
    public:
        /** \brief When comparing doubles if(a < b), do if(a+COMP_MARGIN < b) to avoid
            double precission issues. */
        static constexpr double COMP_MARGIN = std::numeric_limits<double>::epsilon() * 1e5;

         /** \brief Returns true if t1 is at least epsilon-lower than t2, provides robust comparions
             for doubles. */
        static bool isTimeBetterThan
        (double t1, double t2) {
            return t1 + COMP_MARGIN < t2;
        }

        /** \brief An user-implemented absolute value function for integer values. */
        static unsigned int absUI
        (int a) {
            return (a>0) ? (a) : (-a);
        }
};

#endif /* UTILS_H_ */
