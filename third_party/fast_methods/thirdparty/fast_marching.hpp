
#ifndef FAST_MARCHING_H
#define FAST_MARCHING_H

#include <algorithm>
#include <limits>

#include <Image.H>
#include <Globals.h>
#include <PriorityQueue.h>

//const float EPS = 1e-6f;

using namespace Images;

namespace levelset {


        typedef enum { eAlive=0, eTrial=1, eFar=2, eForbidden=3 } eState;


        template<unsigned DIM, typename Pixel = float, int sign = +1>
        class FastMarching : public BaseImage<DIM, unsigned>
        {
        public:

                typedef          BaseImage<DIM, unsigned>    base;
                typedef          BaseImage<DIM, Pixel>       ImageType;
                typedef          BaseImage<DIM, int>         ImageBucket;
                typedef typename base::Index                 Index;
                typedef          PriorityQueue< Index >      CoordQueue;

        

                FastMarching();
                FastMarching(ImageType *_phi);
                virtual ~FastMarching();
                
                void init(ImageType & _phi);

//TODO operator= ?  constructeur copie ?


                
                unsigned nb_alive_points() const;
                const Index& get_alive_point(unsigned i);
                void add_alive_point(const Index& ind);
                void add_trial_point(const Index& ind);
                void add_forbidden_point(const Index& ind);
                void init_trial_from_alive();
                void clear();



                void run();
                void run(Pixel _limit);



        protected:

                Pixel _get_value(const Index& ind) const;
                virtual Pixel _update_value(const Index& ind) const = 0;
                bool _solve_trinome(const Pixel a, const Pixel b, const Pixel c, Pixel &sol_max, Pixel &sol_min) const;

        private:
                void __update_point(const Index& ind);



                Index *              m_alive_tab;    
                unsigned             m_nb_alive;     
                CoordQueue           m_trial_queue;  
                ImageBucket          m_image_bucket; 
        protected:
                Pixel                m_big;          
                ImageType *          m_phi;          

        };


        template <unsigned DIM, typename Pixel, int sign = +1>
        class Eikonal : public FastMarching<DIM, Pixel, sign> 
        {
        public:

                typedef FastMarching<DIM, Pixel, sign> base;
                typedef typename base::ImageType       ImageType;
                typedef typename base::Index           Index;

                Eikonal();
                Eikonal(ImageType * _phi);
                virtual ~Eikonal() {}

        protected:
        
                void set_metric(ImageType * _metric);
                virtual Pixel _update_value(const Index& ind) const;
                Pixel _solve(const Index& ind, Pixel val[DIM]) const;
                
        private:

                ImageType * m_metric;  
        };


        // Constructors - Destructor - Init


        template<unsigned DIM, typename Pixel, int sign>
                FastMarching<DIM, Pixel, sign>::FastMarching()
                : base(), m_phi(NULL), m_alive_tab(NULL), m_nb_alive(0), m_trial_queue(), m_image_bucket(), m_big(std::numeric_limits<Pixel>::max())
        { }

        template<unsigned DIM, typename Pixel, int sign>
                FastMarching<DIM, Pixel, sign>::FastMarching(ImageType *_phi)
                : base(_phi->shape()), m_phi(_phi), m_nb_alive(0), m_trial_queue(), m_image_bucket(_phi->shape()), m_big(std::numeric_limits<Pixel>::max()) 
        {
                m_alive_tab    = new Index[base::size()];
                base::operator = (eFar);
                m_image_bucket = 0;
        }

        template<unsigned DIM, typename Pixel, int sign>
                FastMarching<DIM, Pixel, sign>::~FastMarching() 
        {
                if (m_alive_tab != NULL) delete[] m_alive_tab;
        }

        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::init(ImageType & _phi) 
        {
                // if not the same shape, resize the datas
                if (_phi.shape() != base::shape()){

                        base::resize(_phi.shape());
                        m_image_bucket.resize(_phi.shape());

                        if (m_alive_tab != NULL) delete[] m_alive_tab;
                        m_alive_tab = new Index[base::size()];
                }
                clear();
                m_phi = &_phi;
        }


        //  Access/Modif the Alive-Trial-Forbidden points
        

         template<unsigned DIM, typename Pixel, int sign>
                unsigned FastMarching<DIM, Pixel, sign>::nb_alive_points() const
        {
                return m_nb_alive;
        }
                
        template<unsigned DIM, typename Pixel, int sign>
                const typename FastMarching<DIM, Pixel, sign>::Index& FastMarching<DIM, Pixel, sign>::get_alive_point(unsigned i)
        {
                return m_alive_tab[i];
        }
        
        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::add_alive_point(const Index& ind)
        {
                if (base::operator()(ind) == eFar) {
                        base::operator()(ind) = eAlive;
                        m_alive_tab[m_nb_alive++] = ind;
                }
        }

        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::add_trial_point(const Index& ind)
        {
                // note : the values are not initialized, and they may keep the same value 
                //        (at least for the first one when the fast marching loop begin)

                if (base::operator()(ind) == eFar) {
                        base::operator()(ind) = eTrial;
                        m_image_bucket(ind) = m_trial_queue.push(ind, sign*(*m_phi)(ind));
                }
        }

        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::add_forbidden_point(const Index& ind) 
        {
                if (base::operator()(ind) == eFar) {
                        base::operator()(ind) = eForbidden;
                }
        }

        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::init_trial_from_alive()
        {
                // note : we initialize the values of the points, before inserting them in the queue.
                //        (previous version didn't initialize them)

                // we iterate on the list of alive points, and we add their neighborings point to the trial point heap.


                for (unsigned int j = 0 ; j < nb_alive_points() ; j++)
                {

                        const Index& ind = get_alive_point(j);
                                
                        //std::cout << "alivepoint: " << ind << std::endl;
                        //m_trial_queue.print();

                        for (int i=1 ; i<=DIM ; i++){
                                Index inf = ind;
                                Index sup = ind;
                                inf(i)--;
                                sup(i)++;

                                //std::cout << " voisins " << sup << std::endl;
                                //std::cout << " voisins " << inf << std::endl;

                                if (sup(i)<m_phi->extent(i-1)) __update_point(sup); // add_trial_point(sup);
                                if (inf(i)>=0)                 __update_point(inf); // add_trial_point(inf);
                        }
                }
        }

        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::clear() 
        {
                m_nb_alive = 0;
                base::operator=(eFar);
                m_trial_queue.clear();
                m_image_bucket = 0;
        }


        //  Launch the fast marching


        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::run() 
        {
                run( (sign>0) ? m_big : -m_big );
        }

        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::run(Pixel _limit) 
        {

        //int i = 0;

                // while there are still trial points in the heap
                while(!m_trial_queue.empty()) {

                        // m_trial_queue.print();
                
                        // we get the trial point with the weakest value
                        Index ind = m_trial_queue.top();

                        //std::cout << i++ << std::endl;
                        //std::cout << "top " << ind << ", val = " << sign*(*m_phi)(ind) << std::endl;

                        // If we are at the limit
                        if (sign*(*m_phi)(ind)>=sign*_limit) {
                                // we clear the heap, and mark all its point at 'Far'
                                while (!m_trial_queue.empty()) {
                                        ind = m_trial_queue.top();
                                        m_trial_queue.pop();
                                        base::operator()(ind) = eFar;
                                        m_image_bucket(ind)   = 0;
                                        (*m_phi)(ind)         = _limit;
                                }
                                return;
                        }

                        //std::cout << "pop" << ind << (*m_phi)(ind) << std::endl;

                        // we remove this point from the heap, and add it to the alive points list
                        m_trial_queue.pop();    
                        m_image_bucket(ind) = 0;
                        base::operator()(ind) = eAlive;
                        m_alive_tab[m_nb_alive++] = ind;
                        
                        // we look for the neighborings point of this point, and we update their values
                        for (int i=1 ; i<=DIM ; i++){
                                Index inf = ind;
                                Index sup = ind;
                                inf(i)--;
                                sup(i)++;
                                if (sup(i)<m_phi->extent(i-1)) __update_point(sup);
                                if (inf(i)>=0)                 __update_point(inf);
                        }
                }
        }


        //  member functions


        template<unsigned DIM, typename Pixel, int sign>
                Pixel FastMarching<DIM, Pixel, sign>::_get_value(const Index& ind) const 
        {
                if (base::operator()(ind) <= eTrial)
                        return sign*(*m_phi)(ind);
                return m_big;
        }

        template<unsigned DIM, typename Pixel, int sign>
                bool FastMarching<DIM, Pixel, sign>::_solve_trinome(const Pixel a, const Pixel b, const Pixel c, Pixel &sol_max, Pixel &sol_min) const 
        {
                const Pixel delta = b*b - 4*a*c;
                if (delta < 0) return false;
                const Pixel sqrtDelta = std::sqrt(delta);
                sol_max = (- b + sqrtDelta) / a / 2;
                sol_min = (- b - sqrtDelta) / a / 2;
                return true;
        }

        template<unsigned DIM, typename Pixel, int sign>
                void FastMarching<DIM, Pixel, sign>::__update_point(const Index& ind) 
        {
                //std::cout << "__update_point " << ind << ", val = " << sign*(*m_phi)(ind) << std::endl;

                const eState st = static_cast<eState>(base::operator()(ind));
                if (st == eFar) {
                        // if it was a 'Far' point, we set it to 'Trial', we compute its value and we add it to the heap
                        base::operator()(ind) = eTrial;
                        const Pixel val = _update_value(ind);
                        (*m_phi)(ind) = sign*val;
                //              std::cout << " -----> push" << ind << ", val = " << sign*(*m_phi)(ind) << std::endl;
                        m_image_bucket(ind) = m_trial_queue.push(ind, val);
                }
                else if (st == eTrial) {
                        // if it was a 'Trial' point, we compute its new value
                        const Pixel val = _update_value(ind);
                        //      std::cout << "-----> increase" << ind << ", val = " << sign*(*m_phi)(ind)<< std::endl;
                        // if its value is smaller, we accept it
                        if (val<sign*(*m_phi)(ind)) {
                        //              std::cout << "-----> accept, val = " << val << std::endl;
                                // specify the element in the queue, bucket in the priority queue, and new priority.
                                m_image_bucket(ind) = m_trial_queue.increase_priority(ind, m_image_bucket(ind), val);
                                (*m_phi)(ind) = sign*val;
                        }
                }
        }
        

        // Eikonal Equation


        template<unsigned DIM, typename Pixel, int sign>
                Eikonal<DIM, Pixel, sign>::Eikonal()
                : base(), m_metric(NULL)
        { }

        template<unsigned DIM, typename Pixel, int sign>
                Eikonal<DIM, Pixel, sign>::Eikonal(ImageType *_phi)
                : base(_phi), m_metric(NULL)
        { }
        
        template<unsigned DIM, typename Pixel, int sign>
                void Eikonal<DIM, Pixel, sign>::set_metric(ImageType * _metric)
        {
                m_metric = _metric;
        }

        template<unsigned DIM, typename Pixel, int sign>
                Pixel Eikonal<DIM, Pixel, sign>::_update_value(const Index& ind) const 
        {
                Pixel val[DIM];
                for (int i=1 ; i<=DIM ; i++){
                        Index inf = ind;
                        Index sup = ind;
                        inf(i)--;
                        sup(i)++;

                        //std::cout << ":)" << std::endl;
                        //std::cout << inf << "val=" << (*base::m_phi)(inf) << " fm=" << base::operator()(inf) << std::endl;
                        //std::cout << sup << "val=" << (*base::m_phi)(sup) << " fm=" << base::operator()(sup) << std::endl;
                        
                        if (sup(i)>=base::extent(i-1))  val[i-1] = base::_get_value(inf);
                        else if (inf(i)<0)              val[i-1] = base::_get_value(sup);
                        else val[i-1] = std::min(base::_get_value(inf), base::_get_value(sup));
                }

                //std::cout << "_update_value" << ind ;
                //for (int i = 0 ; i < DIM; i++) std::cout << " " << val[i];
                //std::cout  << std::endl;      
                
                return _solve(ind,val);
        }

        template<unsigned DIM, typename Pixel, int sign>
                Pixel Eikonal<DIM, Pixel, sign>::_solve(const Index& ind, Pixel val[DIM]) const 
        {
                // we reorder value to have val[DIM-1] max
                std::sort(val, val+DIM);

                //std::cout << "_solve" << ind ;
                //for (int i = 0 ; i < DIM; i++) std::cout << " " << val[i];
                //std::cout  << std::endl;
                                
                
                // we initialise the datas of the trinom to solve
                Pixel a, b, c, sol_max, sol_min;
                a = DIM;
                b = 0 ;
                c = (m_metric == NULL)?-1:-Maths::sqr((*m_metric)(ind));
                for (int i=0; i<DIM ; i++){
                         b -= 2*val[i];
                         c += val[i]*val[i];
                }

                //std::cout  << "a=" << a << " b=" << b << " c=" << c << std::endl;

                // we test the differents cases
                for (int i=DIM-1 ; i>0 ; i--){

                        //std::cout  << "DIM=" << i << " a=" << a << " b=" << b << " c=" << c << std::endl;
                        if (val[i]<base::m_big && _solve_trinome(a, b, c, sol_max, sol_min) && sol_max+EPS>=val[i]) {
                                //std::cout << "OK!! solmin" << sol_min <<  " solmax=" << sol_max << std::endl;
                                return sol_max;
                        }
                        //std::cout << "solmin" << sol_min <<  " solmax=" << sol_max << std::endl;
                        a--;
                        b += 2*val[i];
                        c -= val[i]*val[i];
                }
                //std::cout << "last" << val[0] + 1 << std::endl;

                return val[0] + ((m_metric == NULL)? 1 : (*m_metric)(ind));
        }
}


#endif // FAST_MARCHING_H
