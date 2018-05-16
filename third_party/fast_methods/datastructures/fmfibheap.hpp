/*! \class FMDaryHeap
    \brief Wrap for the Boost Fibonacci Heap class to be used in the FM 
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

#ifndef FMFIBHEAP_H_
#define FMFIBHEAP_H_

#include <boost/heap/fibonacci_heap.hpp>

#include <fast_methods/datastructures/fmcompare.hpp>

/// \note for memory efficiency, use map instead of vector for handles_.
template <class cell_t = FMCell> class FMFibHeap {

    /** \brief Shorthand for heap type. */
    typedef boost::heap::fibonacci_heap<const cell_t *, boost::heap::compare<FMCompare<cell_t> > > fib_heap_t;

    /** \brief Shorthand for heap element handle type. */
    typedef typename fib_heap_t::handle_type handle_t;

    public:
        FMFibHeap () {}

        /** \brief Creates a heap with n maximum elements. */
        FMFibHeap (const size_t & n) { handles_.resize(n);}
        virtual ~ FMFibHeap() {}

        /** \brief Sets the maximum number of cells the heap will contain. */
        void setMaxSize
        (const size_t & n) {
            handles_.resize(n);
        }

        /** \brief Pushes a new element into the heap. */
        void push
        (const cell_t * c) {
            handles_[c->getIndex()] = heap_.push(c);
        }

        /** \brief Pops index of the element with lowest value and removes it from the heap. */
        unsigned int popMinIdx
        () {
            const unsigned int idx = heap_.top()->getIndex();
            heap_.pop();
            return idx;
        }

        /** \brief Returns current size of the heap. */
        size_t size
        () const {
            return heap_.size();
        }

        /** \brief Updates the position of the cell in the heap. Its priority can increase or decrease. */
        void update
        (const cell_t * c) {
            heap_.update(handles_[c->getIndex()], c);
        }

        /** \brief Updates the position of the cell in the heap. Its priority can only increase.
            It is more efficient than the update() function if it is ensured that the priority
            will increase. */
        void increase
        (const cell_t * c) {
            heap_.increase(handles_[c->getIndex()],c);
        }
        
        /** \brief Deallocates heap memory. */
        void clear
        () {
            heap_.clear();
            handles_.clear();
        }

        /** \brief Returns true if the heap is empty. */
        bool empty
        () const {
            return heap_.empty();
        }

    protected:
        /** \brief The actual heap for cell_t. */
        fib_heap_t heap_;  /*!< The actual heap for cell_t. */
        
        /** \brief Stores the handles of each cell by keeping the indices: handles_(0) is the handle for
            the cell with index 0 in the grid. Makes possible to update the heap.*/
        std::vector<handle_t> handles_;
};


#endif /* FMFIBHEAP_H_ */

