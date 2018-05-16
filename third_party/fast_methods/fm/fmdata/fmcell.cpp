#include "fmcell.h"
#include "../../console/console.h"

using namespace std;

ostream& operator << 
(ostream & os, const FMCell & c) {
    os << console::str_info("Fast Marching cell information:");
    os << "\t" << "Index: " << c.index_ << '\n'
       << "\t" << "Value: " << c.value_ << '\n'
       << "\t" << "Velocity: " << c.occupancy_ << '\n'
       << "\t" << "State: " ;

    switch (c.state_) {
        case FMState::OPEN:
            os << "OPEN";
            break;
        case FMState::NARROW:
            os << "NARROW";
            break;
        case FMState::FROZEN:
            os << "FROZEN";
            break;
        }
    os << '\n';
    return os;
}

void FMCell::setDefault
() {
    Cell::setDefault();
    value_ = std::numeric_limits<double>::infinity();
    bucket_ = 0;
    hValue_ = 0;
    state_ = FMState::OPEN;
}
