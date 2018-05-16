/*! \class FMUntidyQueue
    \brief Wraps the UntidyQueue implementation by 
    [Jerome Piovano](ftp://ftp-sop.inria.fr/athena/Team/Jerome.Piovano/Doxygen/classlevelset_1_1PriorityQueue.html)
    
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

#ifndef FMUNTIDYQUEUE_HPP_
#define FMUNTIDYQUEUE_HPP_

#include "../../thirdparty/untidy_queue.hpp"
#include "fmcell.h"

/// \todo save buckets here as a hash table instead of saving them in FMCell.
template<class cell_t = FMCell> class FMUntidyQueue {

    public:
        /** \brief Creates an object with s buckets of size s. */
        FMUntidyQueue
        (unsigned s = 1000, double inc = 2) {
            queue_ = new levelset::PriorityQueue<const cell_t * >(s, inc);
        }

        virtual ~FMUntidyQueue() { delete queue_; }

        /** \brief Pushes a new element into the heap. */
        void push
        (cell_t * c) {
            c->setBucket( queue_->push(c, c->getArrivalTime()) );
        }

        /** \brief Returns current size of the heap. */
        size_t size
        () const {
            return queue_->size();
        }

        /** \brief Updates the position of the cell in the priority queue. Its priority can only increase.
             Also updates the bucket of the cell. */
        void increase
        (cell_t * c) {
            c->setBucket( queue_->increase_priority(c, c->getBucket(), c->getArrivalTime()) );
        }

        /** \brief Returns index of the element with \e lowest value (to be popped next). */
        unsigned int topIdx
        (){
            return queue_->top()->getIndex();
        }

        /** \brief Removes the top value of the heap. */
        void pop
        () {
            queue_->pop();
        }

        /** \brief Deallocates heap memory. */
        void clear
        () {
            queue_->clear();
        }

        /** \brief Returns true if the heap is empty. */
        bool empty
        () const {
            return queue_->empty();
        }

    protected:
        /** \brief The actual Unitidy queue for cell_t. */
        levelset::PriorityQueue<const cell_t * > * queue_;
};

#endif /* FMUNTIDYQUEUE_H_ */
