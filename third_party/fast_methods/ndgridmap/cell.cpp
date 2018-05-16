#include "cell.h"
#include "../console/console.h"

using namespace std;

ostream& operator <<
(ostream & os, Cell & c) {
    os << console::str_info("Basic cell information:");
    os << "\t" << "Index: " << c.index_ << '\n'
       << "\t" << "Value: " << c.value_ << '\n'
       << "\t" << "Occupancy: " << c.occupancy_ <<'\n';
    return os;
}

void Cell::setDefault
() {
    value_ = -1;
}
