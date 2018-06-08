#include "hybrid_a_star.h"
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace sdf_tools;

void kinoGridPathFinder::setParameter(double t_ptop_, double t_delta_, double u_max_, double u_delta_, double w_h_, double w_t_, double v_max_, int termination_grid_num_ )
{
    t_prop  = t_ptop_;
    t_delta = t_delta_;
    u_max   = u_max_;
    u_delta = u_delta_;
    w_h     = w_h_;
    w_t     = w_t_;
    v_max   = v_max_;
    termination_grid_num = termination_grid_num_;

    t_steps.clear();
    for(double t = t_delta; t <= t_prop; t += t_delta)
        t_steps.push_back(t);
    
    if( t_steps.back() < t_prop - 0.001 )
        t_steps.push_back(t_prop);
/*
    for(auto ptr:t_steps)
        cout<<ptr<<endl;*/
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
    /*cout<<"check local  map size: "<<X_SIZE<<","<<Y_SIZE<<","<<Z_SIZE<<endl;
    cout<<"check global map size: "<<GLX_SIZE<<","<<GLY_SIZE<<","<<GLZ_SIZE<<endl;
    cout<<"check global origin: "<<gl_xl<<","<<gl_yl<<","<<gl_zl<<endl;*/
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
    closeNodesSequence.clear();
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
    /*double optimal_cost = OptimalControl(node1->state, node2->state).second;
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

    return (1 + tie_breaker) * cost;

//    return (node1->state.head(3)-node2->state.head(3)).norm();

/*    const Vecf<3> dp = node2->state.head(3) - node1->state.head(3);
    const Vecf<3> v0 = node1->state.tail(3);

    decimal_t c1 = -9*dp.dot(dp);
    decimal_t c2 = 12*v0.dot(dp);
    decimal_t c3 = -3*v0.dot(v0);
    decimal_t c4 = 0;
    decimal_t c5 = w_t;

    std::vector<decimal_t> ts = quartic(c5, c4, c3, c2, c1);
    decimal_t t_bar = (node2->state.head(3) - node1->state.head(3)).template lpNorm<Eigen::Infinity>() / v_max;
    ts.push_back(t_bar);

    decimal_t cost = std::numeric_limits<decimal_t>::max();
    for(auto t: ts) {
      if(t < t_bar)
        continue;
      decimal_t c = -c1/3/t/t/t-c2/2/t/t-c3/t + w_t*t;
      if(c < cost)
        cost = c;
    }
    return (1 + tie_breaker) * cost;*/
}

vector<KinoGridNodePtr> kinoGridPathFinder::retrievePath(KinoGridNodePtr current)
{   
    vector<KinoGridNodePtr> path;
    path.push_back(current);

    int cnt = 0;
    while(current->cameFrom != NULL)
    {   
        if(cnt >= 1000)
        {
            ROS_WARN("looping"); ROS_BREAK();
        }
        cout<<current->state(0)<<","<<current->state(1)<<","<<current->state(2)<<endl;
        current = current -> cameFrom;
        path.push_back(current);
        cnt ++;
    }

    ROS_WARN("[hybridAstarSearch] path retrieved");
    return path;
}

vector<Vector3d> kinoGridPathFinder::getClosedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++)
            {   
                if( KinoGridNodeMap[i][j][k]->id == -1 )
                    //visited_nodes.push_back(KinoGridNodeMap[i][j][k]->state.head(3));
                    visited_nodes.push_back( gridIndex2coord(KinoGridNodeMap[i][j][k]->index) );
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

vector<Vector3d> kinoGridPathFinder::getOpenNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++)
            {   
                if( KinoGridNodeMap[i][j][k]->id == 1 )
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
    startPtr -> state.tail(3) = start_vel;
/*    cout<<"startPtr state: \n"<<startPtr->state<<endl;
    cout<<"endPtr state: \n"<<endPtr->state<<endl;*/

    //cout<<"start ptr index: \n"<<startPtr->index<<endl;

    startPtr -> fScore = getHeu(startPtr, endPtr) * w_h;
    openSet.insert( make_pair(startPtr -> fScore, startPtr) ); //put start in open set
    double tentative_gScore;

    num_iter = 0;

    time_in_forward  = 0.0;
    num_ope = 0;

    
    while ( !openSet.empty() && num_iter <= 30000 )
    {   
        num_iter ++;
        currentPtr = openSet.begin() -> second;
        expandedNodes.push_back(currentPtr);
        closeNodesSequence.push_back(gridIndex2coord(currentPtr->index));

        if(    abs(currentPtr->index(0) - endPtr->index(0)) <= termination_grid_num
            && abs(currentPtr->index(1) - endPtr->index(1)) <= termination_grid_num 
            && abs(currentPtr->index(2) - endPtr->index(2)) <= termination_grid_num )
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
            has_path = true;
            //delete neighborPtr; delete currentPtr; delete expandPtr;
            return;
        }         
        openSet.erase(openSet.begin());
        currentPtr -> id = -1; //move current node from open set to closed set.

        Vector3d u;
        VectorXd x0 = currentPtr->state;
        double current_gScore = currentPtr->gScore;
        double s_tie_breaker = 0.0;
        //cout<<"iter:"<<num_iter<<", current state: \n"<<currentPtr -> state<<endl;

        for(double u_x = -u_max; u_x <= u_max; u_x += u_delta )
            for(double u_y = -u_max; u_y <= u_max; u_y += u_delta )
                for(double u_z = -u_max; u_z <= u_max; u_z += u_max )
                {   
                    u << u_x, u_y, u_z;
                    
                    if( forwardSimulation(x0, u, expandPtr) == false) // collision occurs in this transition
                        continue;                    
                    
                    neighborPtr = KinoGridNodeMap[expandPtr->index(0)][expandPtr->index(1)][expandPtr->index(2)];
                    /*if(num_iter < 2)
                    {   
                        ROS_WARN("input");
                        cout<<u<<endl;
                        cout<<"new state: \n"<<expandPtr  -> state<<endl;
                        cout<<"current state: \n"<<currentPtr -> state<<endl;

                        cout<<neighborPtr->index<<endl;
                        cout<<currentPtr->index<<endl;
                    }*/
                    if( neighborPtr->index == currentPtr->index ) //it is still in the same node 
                    {   
                        //cout<<"index :\n"<<neighborPtr->index<<endl;
                        //ROS_BREAK();
                        double tentative_fScore = current_gScore + expandPtr->edge_cost + getHeu(expandPtr, endPtr) * w_h;
                        double real_current_fScore = currentPtr->gScore + getHeu(currentPtr, endPtr) * w_h;

                        if(num_iter < 2)
                        {   
                            cout<<u(0)<<","<<u(1)<<","<<u(2)<<endl;
                            cout<<"current_gScore: "<<current_gScore<<endl;
                            cout<<"expandPtr->edge_cost: "<<expandPtr->edge_cost<<endl;
                            cout<<"Heu of expandPtr :"<<getHeu(expandPtr, endPtr)* w_h<<endl;
                            cout<<"tentative_fScore: "<<tentative_fScore<<endl;
                            cout<<"real_current_fScore fScore: "<<real_current_fScore + s_tie_breaker<<endl;
                        }

                        /*if( tentative_fScore > real_current_fScore + s_tie_breaker)
                            continue;*/
                        if(tentative_fScore <= real_current_fScore + s_tie_breaker)
                        {   
                            //neighborPtr -> input = u;
                            neighborPtr -> state  = expandPtr -> state;
                            neighborPtr -> gScore = current_gScore + expandPtr->edge_cost;
                            neighborPtr -> fScore = current_gScore + expandPtr->edge_cost + getHeu(expandPtr, endPtr) * w_h;
                            
                            if(neighborPtr -> id == 1)
                            {   
                                if(num_iter < 2)
                                    ROS_WARN(" re-insert, take place");
                                
                                openSet.erase(neighborPtr -> nodeMapIt);
                            }

                            neighborPtr -> id = 1;
                            //ROS_WARN(" re-insert");
                            neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                            s_tie_breaker = 0.0;
                        }
                    }
                    else if( neighborPtr -> id != -1 ) //not in closed set 
                    {   
                        tentative_gScore = currentPtr -> gScore + expandPtr->edge_cost; 
                        if(neighborPtr -> id != 1)
                        { //discover a new node
                            neighborPtr -> input = u;
                            neighborPtr -> state = expandPtr -> state;
                            neighborPtr -> id        = 1;
                            neighborPtr -> cameFrom  = currentPtr;
                            neighborPtr -> gScore    = tentative_gScore;
                            neighborPtr -> fScore    = neighborPtr -> gScore + getHeu(neighborPtr, endPtr) * w_h; 
                            neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                        }
                        else if(tentative_gScore < neighborPtr-> gScore)
                        { //in open set and need update
                            neighborPtr -> input = u;
                            neighborPtr -> state = expandPtr -> state;
                            neighborPtr -> cameFrom = currentPtr;
                            neighborPtr -> gScore = tentative_gScore;
                            neighborPtr -> fScore = tentative_gScore + getHeu(neighborPtr, endPtr) * w_h; 
                            openSet.erase(neighborPtr -> nodeMapIt);
                            neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                        }
                    }
                    else
                        continue;
                }
    }

    ros::Time time_2 = ros::Time::now();
    ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
    ROS_WARN("operation num is %d", num_ope );
    has_path = false;
    //delete neighborPtr; delete currentPtr; delete expandPtr;
}

vector<Vector3d> kinoGridPathFinder::getKinoTraj(double resolution)
{     
    vector<Vector3d> statesList;

    KinoGridNodePtr ptr = terminatePtr;
    if( ptr == NULL) // no path found
        return statesList;
    
    VectorXd xt;

    //ROS_WARN("[hybridAstarSearch] check point's index");
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

    ROS_WARN("path size is %d", (int)gridPath.size() );
    for(auto ptr: gridPath)
    {
        path.push_back(ptr->state.head(3));
        /*cout<<"coord: "<<ptr->state(0)<<","<<ptr->state(1)<<","<<ptr->state(2)<<","<<ptr->state(3)<<","<<ptr->state(4)<<","<<ptr->state(5)<<endl;
        cout<<"index: "<<ptr->index(0)<<","<<ptr->index(1)<<","<<ptr->index(2)<<endl;
        cout<<"input: "<<ptr->input(0)<<","<<ptr->input(1)<<","<<ptr->input(2)<<endl;
        cout<<"cost :"<<ptr->gScore<<endl;*/
    }

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

inline bool kinoGridPathFinder::forwardSimulation(VectorXd x0, Vector3d u, KinoGridNodePtr succPtr)
{   
    auto start = std::chrono::high_resolution_clock::now();
    VectorXd xt;
    Vector3i currIdx = coord2gridIndex(x0.head(3));
    int id_x, id_y,id_z;

    for( auto t: t_steps)
    {   
        num_ope ++;
        xt = x0;

        //xt = currentPtr->state;
        getState(xt, u, t);

        if( xt(2) < gl_zl || xt(2) >= gl_zu || xt(0) < gl_xl || xt(0) >= gl_xu || xt(1) < gl_yl || xt(1) >= gl_yu )
            return false;

        coord2gridIndexFast(xt(0), xt(1), xt(2), id_x, id_y, id_z);
        if( KinoGridNodeMap[id_x][id_y][id_z]->occupancy > 0.5) // collision
            return false;
        else if( xt(3) > v_max || xt(4) > v_max || xt(5) > v_max) // velocity too high
            return false;
    }
    double t_exp = t_prop;
    Vector3i succIdx(id_x, id_y, id_z);

   /* if(succIdx == currIdx) // if the success node is in the same node as the current 
    {   
        double t = t_prop;
        //while(succIdx == currentPtr->index && t <= t_max) // try to extend the forwaeding time duration
        while(succIdx == currIdx) // try to extend the forwaeding time duration
        {   
            t += t_delta;

            if( t > t_max)
                return false;
            
            xt = x0;
            getState(xt, u, t);

            if( xt(2) < gl_zl || xt(2) >= gl_zu || xt(0) < gl_xl || xt(0) >= gl_xu || xt(1) < gl_yl || xt(1) >= gl_yu )
                return false;

            coord2gridIndexFast(xt(0), xt(1), xt(2), id_x, id_y, id_z);
            succIdx << id_x, id_y, id_z;

            if( KinoGridNodeMap[id_x][id_y][id_z]->occupancy > 0.5) // collision
                return false;
            else if( xt(3) > v_max || xt(4) > v_max || xt(5) > v_max) // velocity too high
                return false;       
        }

        t_exp = t;
    }*/

    succPtr->index = succIdx;
    succPtr->state = xt;
    succPtr->edge_cost = (u(0) * u(0) + u(1) * u(1) + u(2) * u(2) + w_t ) * t_exp;   

    auto finish = std::chrono::high_resolution_clock::now();
    time_in_forward += std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();                    

    return true;
}