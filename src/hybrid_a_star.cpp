#include "hybrid_a_star.h"
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace sdf_tools;

double getControlCost(double t, VectorXd x0, VectorXd xf)
{   
    double x00, x01, x02, x03, x04, x05;
    double x10, x11, x12, x13, x14, x15;

    x00 = x0(0); x01 = x0(1); x02 = x0(2);
    x03 = x0(3); x04 = x0(4); x05 = x0(5);
    
    x10 = xf(0); x11 = xf(1); x12 = xf(2);
    x13 = xf(3); x14 = xf(4); x15 = xf(5);
    
    double t3 = pow(t, 3);
    double t2 = pow(t, 2);

    double c = t + 4*(x04*x04 + x05*x05 + x13*x13 + x14*x14 + x03*x03 + x03*x13+ x15*x15 + x04*x14 + x05*x15)/t
                 + 12*(x00*x03 + x01*x04 + x02*x05 + x00*x13 - x03*x10 + x01*x14 - x04*x11 + x02*x15 - x05*x12 - x10*x13 - x11*x14 - x12*x15)/t2
                 + 12*(x00*x00 - 2*x00*x10 + x01*x01 + x11*x11 + x02*x02 + x12*x12 + x10*x10 - 2*x01*x11 - 2*x02*x12)/t3;
    return c;
}

int SolveQuartic(const double a, const double b, const double c, const double d, const double e, complex<double>* roots) 
{
    const double a_pw2 = a * a;
    const double b_pw2 = b * b;
    const double a_pw3 = a_pw2 * a;
    const double b_pw3 = b_pw2 * b;
    const double a_pw4 = a_pw3 * a;
    const double b_pw4 = b_pw3 * b;

    const double alpha = -3.0 * b_pw2 / (8.0 * a_pw2) + c / a;
    const double beta  = b_pw3 / (8.0 * a_pw3) - b * c / (2.0 * a_pw2) + d / a;
    const double gamma = -3.0 * b_pw4 / (256.0 * a_pw4) + b_pw2 * c / (16.0 * a_pw3) - b * d / (4.0 * a_pw2) + e / a;

    const double alpha_pw2 = alpha * alpha;
    const double alpha_pw3 = alpha_pw2 * alpha;

    const complex<double> P(-alpha_pw2 / 12.0 - gamma, 0);
    const complex<double> Q(-alpha_pw3 / 108.0 + alpha * gamma / 3.0 - pow(beta, 2.0) / 8.0, 0);
    const complex<double> R = -Q / 2.0 + sqrt(pow(Q, 2.0) / 4.0 + pow(P, 3.0) / 27.0);

    const complex<double> U = pow(R, (1.0 / 3.0));
    complex<double> y;

    const double kEpsilon = 1e-8;
    
    if (fabs(U.real()) < kEpsilon) 
        y = -5.0 * alpha / 6.0 - pow(Q, (1.0 / 3.0));
    else 
        y = -5.0 * alpha / 6.0 - P / (3.0 * U) + U;

    const complex<double> w = sqrt(alpha + 2.0 * y);

    roots[0] = -b / (4.0 * a) +
      0.5 * (w + sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
    roots[1] = -b / (4.0 * a) +
      0.5 * (w - sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));
    roots[2] = -b / (4.0 * a) +
      0.5 * (-w + sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
    roots[3] = -b / (4.0 * a) +
      0.5 * (-w - sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));

    return 4;
}

int SolveQuarticReals(const double a, const double b, const double c, const double d, const double e, const double tolerance, double* roots) 
{
    complex<double> complex_roots[4];
    int num_complex_solutions = SolveQuartic(a, b, c, d, e, complex_roots);
    int num_real_solutions = 0;
    for (int i = 0; i < num_complex_solutions; i++) 
        if (abs(complex_roots[i].imag()) < tolerance) 
          roots[num_real_solutions++] = complex_roots[i].real();

    return num_real_solutions;
}

pair<double, double> OptimalControl(VectorXd x0, VectorXd xf)
{   
/*    cout<<"optimal control: "<<endl;
    cout<<"x0: \n"<<x0<<endl;
    cout<<"xf: \n"<<xf<<endl;*/

    double coef[5];
    double x00, x01, x02, x03, x04, x05;
    double x10, x11, x12, x13, x14, x15;

    x00 = x0(0); x01 = x0(1); x02 = x0(2);
    x03 = x0(3); x04 = x0(4); x05 = x0(5);
    
    x10 = xf(0); x11 = xf(1); x12 = xf(2);
    x13 = xf(3); x14 = xf(4); x15 = xf(5);
             
    // FIX here, all equations should be re-calculated                                                                                                                                     
    coef[0] = 1;
    coef[1] = 0;
    coef[2] = -4*(x03*x03 + x03*x13 + x04*x04 + x04*x14 + x05*x05 + x05*x15 + x13*x13 + x14*x14 + x15*x15);
    coef[3] = 24*(x03*x10 - x01*x04 - x02*x05 - x00*x13 - x00*x03 - x01*x14 + x04*x11 - x02*x15 + x05*x12 + x10*x13 + x11*x14 + x12*x15);
    coef[4] = 36*(-x00*x00 + 2*x00*x10 - x01*x01 + 2*x01*x11 - x02*x02 + 2*x02*x12 - x10*x10 - x11*x11 - x12*x12);
    
    for(int i = 0; i < 5; i++)
        coef[i] /= coef[4];

    double roots[4];
    int root_num = SolveQuarticReals(coef[4], coef[3], coef[2], coef[1], coef[0], 0.0001, roots);
    
    double t_star = -1.0;

    for (int i = 0; i < root_num; i++)
        if(roots[i] > 0.0)
            t_star = 1.0 / roots[i];

    double min_cost = getControlCost(t_star, x0, xf);

    return make_pair(t_star, min_cost);
}

void kinoGridPathFinder::setParameter(double t_ptop_, double t_delta_, double u_max_, double u_delta_, double w_h_, double w_t_, double v_max_ )
{
    t_prop  = t_ptop_;
    t_delta = t_delta_;
    u_max   = u_max_;
    u_delta = u_delta_;
    w_h     = w_h_;
    w_t     = w_t_;
    v_max   = v_max_;

    t_steps.clear();
    for(double t = t_delta; t <= t_prop; t += t_delta)
        t_steps.push_back(t);
    
    if( t_steps.back() < t_prop - 0.001 )
        t_steps.push_back(t_prop);
}

void kinoGridPathFinder::initGridNodeMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    KinoGridNodeMap = new KinoGridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++)
    {
       KinoGridNodeMap[i] = new KinoGridNodePtr * [GLY_SIZE];
       for(int j = 0; j < GLY_SIZE; j++)
       {
           KinoGridNodeMap[i][j] = new KinoGridNodePtr [GLZ_SIZE];
           for( int k = 0; k < GLZ_SIZE;k++)
           {
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                KinoGridNodeMap[i][j][k] = new KinoGridNode(tmpIdx, pos);
           }
       }
    }
}

void kinoGridPathFinder::linkLocalMap(CollisionMapGrid * local_map, Vector3d xyz_l)
{    
    Vector3d coord; 
    for(int64_t i = 0; i < X_SIZE; i++)
    {
        for(int64_t j = 0; j < Y_SIZE; j++)
        {
            for(int64_t k = 0; k < Z_SIZE; k++)
            {   
                coord(0) = xyz_l(0) + (double)(i + 0.5) * resolution;
                coord(1) = xyz_l(1) + (double)(j + 0.5) * resolution;
                coord(2) = xyz_l(2) + (double)(k + 0.5) * resolution;

                Vector3i index = coord2gridIndex(coord);

                if( index(0) >= GLX_SIZE || index(1) >= GLY_SIZE || index(2) >= GLZ_SIZE 
                 || index(0) <  0 || index(1) < 0 || index(2) <  0 )
                    continue;

                KinoGridNodePtr ptr = KinoGridNodeMap[index(0)][index(1)][index(2)];
                ptr->id = 0;
                ptr->occupancy = local_map->Get(i, j, k ).first.occupancy;
            }
        }
    }
}

void kinoGridPathFinder::resetLocalMap()
{   
    //ROS_WARN("expandedNodes size : %d", expandedNodes.size());
    for(auto tmpPtr:expandedNodes)
    {
        tmpPtr->occupancy = 0; // forget the occupancy
        tmpPtr->id = 0;
        tmpPtr->cameFrom = NULL;
        tmpPtr->gScore = inf;
        tmpPtr->fScore = inf;
    }

    for(auto ptr:openSet)
    {   
        KinoGridNodePtr tmpPtr = ptr.second;
        tmpPtr->occupancy = 0; // forget the occupancy
        tmpPtr->id = 0;
        tmpPtr->cameFrom = NULL;
        tmpPtr->gScore = inf;
        tmpPtr->fScore = inf;
    }

    expandedNodes.clear();
    //ROS_WARN("local map reset finish");
}

inline KinoGridNodePtr kinoGridPathFinder::pos2KinoGridNodePtr(Vector3d pos)
{
    Vector3i idx = coord2gridIndex(pos);
    KinoGridNodePtr kino_grid_ptr = new KinoGridNode(idx, pos);

    return kino_grid_ptr;
}

inline Vector3d kinoGridPathFinder::gridIndex2coord(Vector3i index)
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

inline Vector3i kinoGridPathFinder::coord2gridIndex(Vector3d pt)
{
    Vector3i idx;

    //auto start = std::chrono::high_resolution_clock::now();
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);      

    /*auto finish = std::chrono::high_resolution_clock::now();
    time_in_indexing += std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();                    */
    return idx;
}


typedef double decimal_t;
template <int N> 
using Vecf = Eigen::Matrix<decimal_t, N, 1>;

/* **************************************************************** */
inline std::vector<decimal_t> cubic(decimal_t a, decimal_t b, decimal_t c, decimal_t d) {
  std::vector<decimal_t> dts;

  decimal_t a2 = b / a;
  decimal_t a1 = c / a;
  decimal_t a0 = d / a;
  //printf("a: %f, b: %f, c: %f, d: %f\n", a, b, c, d);

  decimal_t Q = (3*a1-a2*a2)/9;
  decimal_t R = (9*a1*a2-27*a0-2*a2*a2*a2)/54;
  decimal_t D = Q*Q*Q + R*R;
  //printf("R: %f, Q: %f, D: %f\n", R, Q, D);
  if(D > 0) {
    decimal_t S = std::cbrt(R+sqrt(D));
    decimal_t T = std::cbrt(R-sqrt(D));
    //printf("S: %f, T: %f\n", S, T);
    dts.push_back(-a2/3+(S+T));
    return dts;
  }
  else if(D == 0) {
    decimal_t S = std::cbrt(R);
    dts.push_back(-a2/3+S+S);
    dts.push_back(-a2/3-S);
    return dts;
  }
  else {
    decimal_t theta = acos(R/sqrt(-Q*Q*Q));
    dts.push_back(2*sqrt(-Q)*cos(theta/3)-a2/3);
    dts.push_back(2*sqrt(-Q)*cos((theta+2*M_PI)/3)-a2/3);
    dts.push_back(2*sqrt(-Q)*cos((theta+4*M_PI)/3)-a2/3);
    return dts;
  }
}

/* **************************************************************** */
inline std::vector<decimal_t> quartic(decimal_t a, decimal_t b, decimal_t c, decimal_t d, decimal_t e) {
  std::vector<decimal_t> dts;

  decimal_t a3 = b / a;
  decimal_t a2 = c / a;
  decimal_t a1 = d / a;
  decimal_t a0 = e / a;

  std::vector<decimal_t> ys = cubic(1, -a2, a1*a3-4*a0, 4*a2*a0-a1*a1-a3*a3*a0);
  decimal_t y1 = ys.front();
  //printf("y1: %f\n", y1);
  decimal_t r = a3*a3/4-a2+y1;
  //printf("r: %f\n", r);

  //printf("a = %f, b = %f, c = %f, d = %f, e = %f\n", a, b, c, d, e);
  if(r < 0)
    return dts;

  decimal_t R = sqrt(r);
  decimal_t D, E;
  if(R != 0) {
    D = sqrt(0.75*a3*a3-R*R-2*a2+0.25*(4*a3*a2-8*a1-a3*a3*a3)/R);
    E = sqrt(0.75*a3*a3-R*R-2*a2-0.25*(4*a3*a2-8*a1-a3*a3*a3)/R);
  }
  else {
    D = sqrt(0.75*a3*a3-2*a2+2*sqrt(y1*y1-4*a0));
    E = sqrt(0.75*a3*a3-2*a2-2*sqrt(y1*y1-4*a0));
  }

  if(!std::isnan(D)) {
    dts.push_back(-a3/4+R/2+D/2);
    dts.push_back(-a3/4+R/2-D/2);
  }
  if(!std::isnan(E)) {
    dts.push_back(-a3/4-R/2+E/2);
    dts.push_back(-a3/4-R/2-E/2);
  }

  return dts;
}

inline double kinoGridPathFinder::getHeu(KinoGridNodePtr node1, KinoGridNodePtr node2)
{   
/*    double optimal_cost = OptimalControl(node1->state, node2->state).second;
    return tie_breaker * optimal_cost;*/

    const Vecf<3> dp = node2->state.head(3) - node1->state.head(3);
    const Vecf<3> v0 = node1->state.tail(3);
    const Vecf<3> v1 = node2->state.tail(3);

    decimal_t c1 = -36*dp.dot(dp);
    decimal_t c2 = 24*(v0+v1).dot(dp);
    decimal_t c3 = -4*(v0.dot(v0)+v0.dot(v1)+v1.dot(v1));
    decimal_t c4 = 0;
    decimal_t c5 = w_t;

    std::vector<decimal_t> ts = quartic(c5, c4, c3, c2, c1);
    decimal_t t_bar = (node1->state.head(3) - node2->state.head(3)).template lpNorm<Eigen::Infinity>() / v_max;
    ts.push_back(t_bar);

    decimal_t cost = std::numeric_limits<decimal_t>::max();
    for(auto t: ts) 
    {
      if(t < t_bar)
        continue;
      decimal_t c = -c1/3/t/t/t-c2/2/t/t-c3/t + w_t * t;
      if(c < cost)
        cost = c;
    }

/*    const Vecf<3> dp = node2->state.head(3) - node1->state.head(3);
    const Vecf<3> v0 = node1->state.tail(3);

    decimal_t c1 = -9*dp.dot(dp);
    decimal_t c2 = 12*v0.dot(dp);
    decimal_t c3 = -3*v0.dot(v0);
    decimal_t c4 = 0;
    decimal_t c5 = 1.0;

    std::vector<decimal_t> ts = quartic(c5, c4, c3, c2, c1);
    decimal_t t_bar = (node2->state.head(3) - node1->state.head(3)).template lpNorm<Eigen::Infinity>() / v_max;
    ts.push_back(t_bar);

    decimal_t cost = std::numeric_limits<decimal_t>::max();
    for(auto t: ts) {
      if(t < t_bar)
        continue;
      decimal_t c = -c1/3/t/t/t-c2/2/t/t-c3/t + 1.0*t;
      if(c < cost)
        cost = c;
    }
*/
    return tie_breaker * cost * w_h;
}

vector<KinoGridNodePtr> kinoGridPathFinder::retrievePath(KinoGridNodePtr current)
{   
    vector<KinoGridNodePtr> path;
    path.push_back(current);

    while(current->cameFrom != NULL)
    {   
        //cout<<current->state(0)<<","<<current->state(1)<<","<<current->state(2)<<endl;
        current = current -> cameFrom;
        path.push_back(current);
    }

    return path;
}

vector<Vector3d> kinoGridPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++)
            {   
                if( KinoGridNodeMap[i][j][k]->id != 0)
                    //visited_nodes.push_back(KinoGridNodeMap[i][j][k]->state.head(3));
                    visited_nodes.push_back( gridIndex2coord(KinoGridNodeMap[i][j][k]->index) );
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

void kinoGridPathFinder::hybridAstarSearch(Vector3d start_pt, Vector3d start_vel, Vector3d end_pt, Vector3d end_vel)
{   
    ros::Time time_1 = ros::Time::now();    
    
    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx = coord2gridIndex(end_pt);

    KinoGridNodePtr startPtr = KinoGridNodeMap[start_idx(0)][start_idx(1)][start_idx(2)];
    KinoGridNodePtr endPtr   = KinoGridNodeMap[end_idx(0)][end_idx(1)][end_idx(2)];
    endPtr->state.tail(3) = end_vel;

    openSet.clear();

    KinoGridNodePtr neighborPtr  = NULL;
    KinoGridNodePtr currentPtr   = NULL;
    KinoGridNodePtr expandPtr = new KinoGridNode();
    terminatePtr = NULL;

    startPtr -> gScore = 0;
    startPtr -> id = 1; //put start node in open set
    startPtr -> state.head(3) = start_pt;
    startPtr -> state.tail(3) = start_vel;
    startPtr -> fScore = getHeu(startPtr, endPtr);

    openSet.insert( make_pair(startPtr -> fScore, startPtr) ); //put start in open set
    double tentative_gScore;

    int num_iter = 0;

    time_in_forward  = 0.0;
    num_ope = 0;

    while ( !openSet.empty() && num_iter <= 30000 )
    {   
        num_iter ++;
        currentPtr = openSet.begin() -> second;
        if(sqrt( (currentPtr->state(0)-endPtr->state(0)) * (currentPtr->state(0)-endPtr->state(0)) + 
                 (currentPtr->state(1)-endPtr->state(1)) * (currentPtr->state(1)-endPtr->state(1)) + 
                 (currentPtr->state(2)-endPtr->state(2)) * (currentPtr->state(2)-endPtr->state(2)) ) <= 1.0 )
        /*&& sqrt( (current->state(3)-endPtr->state(3)) * (current->state(3)-endPtr->state(3)) + 
                 (current->state(4)-endPtr->state(4)) * (current->state(4)-endPtr->state(4)) + 
                 (current->state(5)-endPtr->state(5)) * (current->state(5)-endPtr->state(5)) ) <= 1.0) */
        {
            ROS_WARN("[Astar]Reach goal..");
            //cout << "goal coord: " << endl << current->real_coord << endl; 
            cout << "total number of iteration used in Astar: " << num_iter  << endl;
            ros::Time time_2 = ros::Time::now();
            ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
            ROS_WARN("Time consume in forward simulation is %f", time_in_forward * 1.0e-9 );
            ROS_WARN("operation num is %d", num_ope );
            gridPath = retrievePath(currentPtr);
            terminatePtr = currentPtr;
            //delete neighborPtr; delete currentPtr; delete expandPtr;
            return;
        }         
        openSet.erase(openSet.begin());
        currentPtr -> id = -1; //move current node from open set to closed set.
        expandedNodes.push_back(currentPtr);

        Vector3d u;
        for(double u_x = -u_max; u_x <= u_max; u_x += u_delta )
            for(double u_y = -u_max; u_y <= u_max; u_y += u_delta )
                for(double u_z = -u_max; u_z <= u_max; u_z += u_max )
                {   
                    u << u_x, u_y, u_z;
                    if( forwardSimulation(currentPtr, u, expandPtr) == false) // collision occurs in this transition
                        continue;                    
                    
                    neighborPtr = KinoGridNodeMap[expandPtr->index(0)][expandPtr->index(1)][expandPtr->index(2)];
                    if( neighborPtr -> id != -1 ) //not in closed set or it is still in the same node 
                    {
                        tentative_gScore = currentPtr -> gScore + expandPtr->edge_cost; 
                        if(neighborPtr -> id != 1)
                        { //discover a new node
                            neighborPtr -> input = u;
                            neighborPtr -> state = expandPtr -> state;
                            neighborPtr -> id        = 1;
                            neighborPtr -> cameFrom  = currentPtr;
                            neighborPtr -> gScore    = tentative_gScore;
                            neighborPtr -> fScore    = neighborPtr -> gScore + getHeu(neighborPtr, endPtr); 
                            neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                        }
                        else if(tentative_gScore < neighborPtr-> gScore)
                        { //in open set and need update
                            neighborPtr -> input = u;
                            neighborPtr -> state = expandPtr -> state;
                            neighborPtr -> cameFrom = currentPtr;
                            neighborPtr -> gScore = tentative_gScore;
                            neighborPtr -> fScore = tentative_gScore + getHeu(neighborPtr, endPtr); 
                            openSet.erase(neighborPtr -> nodeMapIt);
                            neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                        }
                    }
                }
    }

    ros::Time time_2 = ros::Time::now();
    ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
    ROS_WARN("operation num is %d", num_ope );
    //delete neighborPtr; delete currentPtr; delete expandPtr;
}

vector<Vector3d> kinoGridPathFinder::getKinoTraj(double resolution)
{     
    vector<Vector3d> statesList;

    KinoGridNodePtr ptr = terminatePtr;
    if( ptr == NULL) // no path found
        return statesList;
    
    VectorXd xt;

    while( ptr->cameFrom != NULL ) 
    {   
        Vector3d u = ptr->input;    
        for(double t = t_prop; t >= resolution; t -= resolution)
        {
            xt = ptr->cameFrom->state;
            getState( xt, u, t );
            statesList.push_back( xt.head(3) );
        }

        ptr = ptr->cameFrom;
    }

    reverse(statesList.begin(), statesList.end());
    return statesList;
}

vector<Vector3d> kinoGridPathFinder::getPath()
{   
    vector<Vector3d> path;

    for(auto ptr: gridPath)
        path.push_back(ptr->state.head(3));

    reverse(path.begin(), path.end());
    return path;
}

void kinoGridPathFinder::resetPath()
{
    gridPath.clear();
}

inline void kinoGridPathFinder::getState( VectorXd & xt, Vector3d u, double t )
{
    double t2 = t * t / 2.0;

    xt(0) += t * xt(3) + t2 * u(0);
    xt(1) += t * xt(4) + t2 * u(1);
    xt(2) += t * xt(5) + t2 * u(2);

    xt(3) += t * u(0);
    xt(4) += t * u(1);
    xt(5) += t * u(2);
}

inline void kinoGridPathFinder::coord2gridIndexFast(double x, double y, double z, int & id_x, int & id_y, int & id_z)
{
    id_x = int( (x - gl_xl) * inv_resolution);
    id_y = int( (y - gl_yl) * inv_resolution);
    id_z = int( (z - gl_zl) * inv_resolution);      
}

inline bool kinoGridPathFinder::forwardSimulation(KinoGridNodePtr p_cur, Vector3d u, KinoGridNodePtr succPtr)
{   
    auto start = std::chrono::high_resolution_clock::now();
    VectorXd xt;
    int id_x, id_y,id_z;

    for( auto t: t_steps)
    {   
        num_ope ++;

        xt = p_cur->state;
        getState(xt, u, t);

        if( xt(2) < gl_zl || xt(2) >= gl_zu || xt(0) < gl_xl || xt(0) >= gl_xu || xt(1) < gl_yl || xt(1) >= gl_yu )
            return false;

        coord2gridIndexFast(xt(0), xt(1), xt(2), id_x, id_y, id_z);
        if( KinoGridNodeMap[id_x][id_y][id_z]->occupancy > 0.5) // collision
            return false;
        else if( xt(3) > v_max || xt(4) > v_max || xt(5) > v_max) // velocity too high
            return false;

    }

    succPtr->index = Vector3i(id_x, id_y, id_z);
    succPtr->state = xt;
    succPtr->edge_cost = u(0) * u(0) + u(1) * u(1) + u(2) * u(2) + w_t * t_prop;

    auto finish = std::chrono::high_resolution_clock::now();
    time_in_forward += std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();                    

    return true;
}