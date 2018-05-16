 #ifndef PRIORITY_QUEUE_H
#define PRIORITY_QUEUE_H

#include <list>
#include <cmath>
#include <iomanip>

namespace levelset {



//TODO expliquer que c est en O(n) mais au detriment de la stabilite (minime) si le range des valeur contenues simultanement dans la queue est > incmax.




//TODO tester sans debug car sur 2000x2000 1000buckets nous => 11sec eux, 1.5 sec. (a la limite remettre les iters + buckets , pour remove temps constant, et tester complexite en augmentant taille progressivement.)\
//TODO en fait en utilisant gprof, increase priority est appele 4000 fois sur les 4000000 pop et push ==> iterateur negligeable pour le removede plus moyenne remplisage = 13 avec 257 max elet par list donc ...
//TODO 3sec en -O4 sans debugs.6sec en -O4 version Jphil .... 
//TODO tester les increments >> 1.5 (wtf?)

//TODO compris les >> increment : quabnd on fait init trial from alive, recupere les anciennes valeur de phi, du coup fais des push successifs, sans pop, et ptet que 2 push separes de 1
// fortement possible car pas de PDEreinit
// TODO finalemenmt implemente le push et increase pas dans l'ordre . aucun tests possibles sur les increments ,  mais en pratique (pour le fast marching) devrait marcher.
// TODO juste preciser que l ensemble des valeurs de la queue doit etre comprise dans un range precis.

//TODO ATTENTION !!!!!!!!!!!!!!!!!! si precision float , desfois incmax > 1.05 , si precision double passe avec incmax = 1.002 . en theorie <=1 pas plus.
//TODO rendre + clair la partie push dans le desordre.


//TODO incmax chelou a prendre, car meme si a 2, plante pour le cas d utilisation me masque (parties trop pointues.)

//TODO rq : kan laisse passser si bug, marche quand meme au top., meme avec 1 incmax, si foire a 2.

        template <class Element>
        class PriorityQueue
        {
        public:
                typedef std::list<Element>               ElementList;


                PriorityQueue(unsigned _size    = 1000,   
                              double   _inc_max = 2       
                        );

                ~PriorityQueue();



                bool empty() const;
                int push(const Element &e, double t);
                void pop();
                const Element & top();
                int increase_priority(const Element& e, int bucket, double t_new);
                void clear();
                unsigned size() const;

                void print() const;

        private:



                ElementList *    m_tab;      
                unsigned         m_size;     
                unsigned         m_nb_elem;  
                double           m_t0;       
                int              m_i0;       
                double           m_delta;    
                double           m_inc_max;  
        };


        // Constructors - Destructor


        template <class Element>
                PriorityQueue<Element>::PriorityQueue(unsigned _size, double _inc_max) 
                : m_size(_size), m_t0(0), m_i0(0), m_nb_elem(0),m_inc_max(_inc_max)
        {
                m_tab   = new ElementList[_size];
                m_delta = _inc_max / static_cast<double>(_size);
        }

        template <class Element>
                PriorityQueue<Element>::~PriorityQueue() 
        {
                delete[] m_tab; 
        }


        // Access/Modif the queue


        template <class Element>
                bool PriorityQueue<Element>::empty() const
        {
                return (m_nb_elem == 0);
        }
                
        template <class Element>
                int PriorityQueue<Element>::push(const Element &e, double t)
        {
                int i = 0;

                if (empty()) {
                        m_t0    = t;
                        m_i0    = 0;
                }
                else {

                        // test to prevent errors
                        if (t-m_t0 > m_inc_max) {
                                std::cerr << "Error : Increment=" << t-m_t0 << " >> Inc max" << m_inc_max << std::endl;
                                std::cerr << "        e=" << e << " , t=" << t << " , m_t0=" << m_t0 << std::endl;
                                //print();
                                //TODO laisse passer le prog kan mm en cas de bug.
                                //exit(2);
                        }

                        // compute the indice in the circular array.
                        i = static_cast<int>(floor((t-m_t0) * m_size / m_inc_max));

                        // compute the true indice in the regular array.
                        if (i >= 0){
                                i = (i >= m_size)? m_size-1 : i;
                                i = (i + m_i0) % m_size;
                        }
                        else{
                                // if t < m_t0, compute new values for m_i0, m_t0, and verify that the end of the circular array is empty of at least |i| bucket.
                                while (i++ < 0){

                                        m_i0 = (m_i0-1 >= 0)? m_i0-1 : m_i0-1+m_size;
                                        m_t0 -= m_delta;

                                        if (!m_tab[m_i0].empty()) {
                                                std::cerr << "Error : While shifting the begining of the queue, circular array not empty at the end" << std::endl;
                                                std::cerr << "        e=" << e << " , t=" << t << std::endl;
                                                //print();
                                                // TODO laisse passer le prog kan mm en cas de bug, insere le point ds le premier bucket en  premiere position, meme si n y appartien pas vraiment.
                                                m_i0 = (m_i0+1) % m_size;
                                                m_t0 += m_delta;
                                                break;
                                                //exit(2);
                                        }
                                }
                                i = m_i0;
                        }
                }

                m_nb_elem++;
                m_tab[i].push_front(e);

                return i;
        }
                
        template <class Element>
                void PriorityQueue<Element>::pop()
        {
                if (empty()) return;

                // shift m_i0 in order to be on a non empty level of quantization
                while (m_tab[m_i0].empty()){
                        m_i0 = (m_i0 + 1) % m_size;
                        m_t0 += m_delta;
                }

                m_nb_elem--;
                m_tab[m_i0].pop_back();

        }

        template <class Element>
                const Element & PriorityQueue<Element>::top()
        {
                // shift m_i0 in order to be on a non empty level of quantization
                while (m_tab[m_i0].empty()){
                        m_i0 = (m_i0 + 1) % m_size;
                        m_t0 += m_delta;
                }

                return m_tab[m_i0].back();
        }
                
        template <class Element>
                int PriorityQueue<Element>::increase_priority(const Element& e, int bucket, double t_new)
        {
                if (fabs(t_new-m_t0) > m_inc_max) {
                        std::cerr << "Error : Increment=" << t_new-m_t0 << " >> Inc max" << m_inc_max << std::endl;
                        std::cerr << "        e=" << e << " , t=" << t_new << " , m_t0=" << m_t0 << std::endl;
                        //print();
                        //TODO laisse passer le prog kan mm en cas de bug.
                        //exit(2);
                }
                
                // compute the indice in the circular array
                int  i = static_cast<int>(floor((t_new-m_t0) * m_size / m_inc_max));

                // compute the true indice in the regular array.
                if (i >= 0){
                        i = (i >= m_size)? m_size-1 : i;
                        i = (i + m_i0) % m_size;
                }
                else {
                        // if t < m_t0, compute new values for m_i0, m_t0, and verify that the end of the circular array is empty of at least |i| bucket.
                        while (i++ < 0){

                                m_i0 = (m_i0-1 >= 0)? m_i0-1 : m_i0-1+m_size;
                                m_t0 -= m_delta;

                                if (!m_tab[m_i0].empty()) {
                                        std::cerr << "Error : While shifting the begining of the queue, circular array not empty at the end" << std::endl;
                                        std::cerr << "        e=" << e << " , t=" << t_new << std::endl;
                                        //print();
                                        // TODO laisse passer le prog kan mm en cas de bug, insere le point ds le premier bucket en  premiere position, meme si n y appartien pas vraiment.
                                        m_i0 = (m_i0+1) % m_size;
                                        m_t0 += m_delta;
                                        break;
                                        //exit(2);
                                }
                        }
                        i = m_i0;
                }

                m_tab[bucket].remove(e);
                m_tab[i].push_front(e);

                return i;
        }
                
        template <class Element>
                void PriorityQueue<Element>::clear()
        {
                m_nb_elem = 0;
                m_t0      = 0;
                m_i0      = 0;
                for (int i=0 ; i<m_size ; i++)
                        m_tab[i].clear();
        }

        template <class Element>
                unsigned PriorityQueue<Element>::size() const
        {
                return m_nb_elem;
        }

        template <class Element>
                void PriorityQueue<Element>::print() const
        {
                std::cout << "-------------------------------------"  << std::endl;
                std::cout << "m_t0 = " << m_t0 << std::endl;
                std::cout << "m_i0 = " << m_i0 << std::endl;
                std::cout << "size = " << m_nb_elem << std::endl;
                for (int i = 0; i < m_size; i++){
                        if (!m_tab[i].empty()){
                                // index in the circular array, starting at m_i0
                                int j = (i-(int)m_i0 >= 0)? i-(int)m_i0 : i-(int)m_i0+m_size;
        
                                // print infos concerning the bucket
                                std::cout << "  i=" << std::setw(2) << i 
                                          << "| l=" << std::setw(2) << j 
                                          << "| p=" << std::setw(5) << m_t0 + j*m_delta 
                                          << "->"   << std::setw(5) << m_t0 + j*m_delta + m_delta;
                                // print the elements in the bucket
                                for (typename ElementList::iterator it = m_tab[i].begin() ; it!= m_tab[i].end() ; it++)
                                        std::cout << " " << *it ;
        
                                std::cout << std::endl; 
                        }
                }
        }

}

#endif // PRIORITY_QUEUE_H
