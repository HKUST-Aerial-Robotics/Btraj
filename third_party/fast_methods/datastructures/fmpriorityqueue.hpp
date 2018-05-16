/*! \class FMPriorityQueue
    \brief Wrap for the Boost Priority Queue class to be used in the FM 
    algorithms. Ready to be used with FMCell and derived types.
    
    Copyright (C) 2014 Javier V. Gomez and Jose Pardeiro
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

#ifndef FMPRIORITYQUEUE_H_
#define FMPRIORITYQUEUE_H_

#include <boost/heap/priority_queue.hpp>

#include <fast_methods/datastructures/fmcompare.hpp>

template <class cell_t = FMCell> class FMPriorityQueue{

    public:
        FMPriorityQueue () {}

        /** \brief Shorthand for heap element handle type. */
        FMPriorityQueue (const int & n) { heap_.reserve(n); }

        virtual ~ FMPriorityQueue() {}

        /** \brief Sets the maximum number of cells the heap will contain. */
        void setMaxSize
        (const int & n) {
            heap_.reserve(n);
        }

        /** \brief Pushes a new element into the heap. */
        void push 
        (const cell_t * c) {
            heap_.push(c);
        }

        /** \brief Priority queues do not allow key increasing. Therefore, it pushes the element again.
             This is done so that SFMM is implemented as FMM with this heap. */
        void increase
        (const cell_t * c) {
            heap_.push(c);
        }

        /** \brief Pops index of the element with lowest value and removes it from the heap. */ 
        int popMinIdx
        () {
            const int idx = heap_.top()->getIndex();
            heap_.pop();
            return idx;
        }

        /** \brief Returns current size of the heap. */
        size_t size
        () const {
            return heap_.size();
        }

        /** \brief Deallocates heap memory. */
        void clear
        () {
            heap_.clear();
        }

        /** \brief Returns true if the heap is empty. */
        bool empty
        () const {
            return heap_.empty();
        }

    protected:
        /** \brief The actual queue for FMCells. */
        boost::heap::priority_queue<const cell_t *, boost::heap::compare<FMCompare<cell_t> > > heap_;
};


#endif /* FMPRIORITYQUEUE_H_ */
