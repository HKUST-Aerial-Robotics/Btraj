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

inline vector<double> cubic(double a, double b, double c, double d) 
{
    vector<double> dts;

    double a2 = b / a;
    double a1 = c / a;
    double a0 = d / a;

    double Q = (3*a1-a2*a2)/9;
    double R = (9*a1*a2-27*a0-2*a2*a2*a2)/54;
    double D = Q*Q*Q + R*R;
    if(D > 0) 
    {
        double S = std::cbrt(R+sqrt(D));
        double T = std::cbrt(R-sqrt(D));
        dts.push_back(-a2/3+(S+T));
        return dts;
    }
    else if(D == 0) 
    {
        double S = std::cbrt(R);
        dts.push_back(-a2/3+S+S);
        dts.push_back(-a2/3-S);
        return dts;
    }
    else 
    {
        double theta = acos(R/sqrt(-Q*Q*Q));
        dts.push_back(2*sqrt(-Q)*cos(theta/3)-a2/3);
        dts.push_back(2*sqrt(-Q)*cos((theta+2*M_PI)/3)-a2/3);
        dts.push_back(2*sqrt(-Q)*cos((theta+4*M_PI)/3)-a2/3);
        return dts;
    }
}

inline vector<double> quartic(double a, double b, double c, double d, double e) 
{
    vector<double> dts;

    double a3 = b / a;
    double a2 = c / a;
    double a1 = d / a;
    double a0 = e / a;

    vector<double> ys = cubic(1, -a2, a1*a3-4*a0, 4*a2*a0-a1*a1-a3*a3*a0);
    double y1 = ys.front();
    double r = a3*a3/4-a2+y1;
    if(r < 0)
        return dts;

    double R = sqrt(r);
    double D, E;
    if(R != 0) 
    {
        D = sqrt(0.75*a3*a3-R*R-2*a2+0.25*(4*a3*a2-8*a1-a3*a3*a3)/R);
        E = sqrt(0.75*a3*a3-R*R-2*a2-0.25*(4*a3*a2-8*a1-a3*a3*a3)/R);
    }
    else 
    {
        D = sqrt(0.75*a3*a3-2*a2+2*sqrt(y1*y1-4*a0));
        E = sqrt(0.75*a3*a3-2*a2-2*sqrt(y1*y1-4*a0));
    }

    if(!std::isnan(D)) 
    {
        dts.push_back(-a3/4+R/2+D/2);
        dts.push_back(-a3/4+R/2-D/2);
    }
    if(!std::isnan(E)) 
    {
        dts.push_back(-a3/4-R/2+E/2);
        dts.push_back(-a3/4-R/2-E/2);
    }

  return dts;
}

inline double kinoGridPathFinder::getHeu(KinoGridNodePtr node1, KinoGridNodePtr node2)
{   
    const Vector3d dp = node2->state.head(3) - node1->state.head(3);
    const Vector3d v0 = node1->state.tail(3);
    const Vector3d v1 = node2->state.tail(3);

    double c1 = -36*dp.dot(dp);
    double c2 = 24*(v0+v1).dot(dp);
    double c3 = -4*(v0.dot(v0)+v0.dot(v1)+v1.dot(v1));
    double c4 = 0;
    double c5 = w_t;

    std::vector<double> ts = quartic(c5, c4, c3, c2, c1);
    double t_bar = (node1->state.head(3) - node2->state.head(3)).lpNorm<Infinity>() / v_max;
    ts.push_back(t_bar);

    double cost = std::numeric_limits<double>::max();
    double t_d  = t_bar;

    for(auto t: ts) 
    {
        if(t < t_bar)
            continue;
        double c = -c1/3/t/t/t-c2/2/t/t-c3/t + w_t * t;
        if(c < cost)
        {
            cost = c;
            t_d = t;
        }
    }

    node1->optimal_t = t_d;
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

inline bool kinoGridPathFinder::shotHeu(KinoGridNodePtr node1, KinoGridNodePtr node2)
{   
    cnt_shot ++;
    if(cnt_shot < N_max_shot)
    {
        return false;
    }
    else
    {   
        cnt_shot = 0;
        N_max_shot = ceil( (node1->state.head(3) - node2->state.head(3)).norm() / dis_shot * N_max );
        cout<<"N_max_shot: "<<N_max_shot<<endl;
    }

    const Vector3d p0 = node1->state.head(3);
    const Vector3d dp = node2->state.head(3) - p0;
    const Vector3d v0 = node1->state.tail(3);
    const Vector3d v1 = node2->state.tail(3);
    const Vector3d dv = v1 - v0;
    const double t_d = node1->optimal_t;

/*    double c1 = -36*dp.dot(dp);
    double c2 = 24*(v0+v1).dot(dp);
    double c3 = -4*(v0.dot(v0)+v0.dot(v1)+v1.dot(v1));
    double c4 = 0;
    double c5 = w_t;

    std::vector<double> ts = quartic(c5, c4, c3, c2, c1);
    double t_bar = (node1->state.head(3) - node2->state.head(3)).lpNorm<Infinity>() / v_max;
    ts.push_back(t_bar);

    double cost = std::numeric_limits<double>::max();
    double t_d  = t_bar;

    for(auto t: ts) 
    {
        if(t < t_bar)
            continue;
        double c = -c1/3/t/t/t-c2/2/t/t-c3/t + w_t * t;
        if(c < cost)
        {
            cost = c;
            t_d = t;
        }
    }*/

//  ****** now check the feasibility of the optimal polynomial, by using t_d

    Vector3d a = 1.0 / 6.0 * ( -12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv );
    Vector3d b =    0.5 *    (  6.0  / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv );
    MatrixXd coef(3, 4);

    coef.col(3) = a;
    coef.col(2) = b;
    coef.col(1) = v0;
    coef.col(0) = p0;

// *** the OPTIMAL polynomial is : 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0 * t + p0; denote as : a*t^3 + b*t^2 + v0*t + p0
    Vector3d coord;
    VectorXd poly1d, t;
    int id_x, id_y,id_z;
    
    for(double time = t_delta; time <= t_d; time += t_delta )
    {   
        t = VectorXd::Zero(4); 
        for(int j = 0; j < 4; j ++)
            t(j) = pow(time, j);

        for ( int dim = 0; dim < 3; dim++ )
        {
            poly1d = coef.row(dim);
            coord(dim) = poly1d.dot(t);
        }

        if( coord(2) < gl_zl || coord(2) >= gl_zu || coord(0) < gl_xl || coord(0) >= gl_xu || coord(1) < gl_yl || coord(1) >= gl_yu )
            return false; //continue;

        coord2gridIndexFast(coord(0), coord(1), coord(2), id_x, id_y, id_z);
        if( KinoGridNodeMap[id_x][id_y][id_z]->occupancy > 0.5) // collision
            return false;
    }

    coef_shot = coef;
    t_shot    = t_d;
    is_shot_succ = true;
    return true;
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

    startPtr -> fScore = getHeu(startPtr, endPtr) * w_h;
    openSet.insert( make_pair(startPtr -> fScore, startPtr) ); //put start in open set
    double tentative_gScore;

    num_iter = 0;
    num_ope = 0;
    time_in_forward  = 0.0;
    is_shot_succ = false;
    N_max_shot = 10;
    cnt_shot = 0;
    dis_shot = (startPtr->state.head(3) - endPtr->state.head(3)).norm();

    while ( !openSet.empty() && num_iter <= 30000 )
    {   
        num_iter ++;
        currentPtr = openSet.begin() -> second;
        expandedNodes.push_back(currentPtr);
        closeNodesSequence.push_back(gridIndex2coord(currentPtr->index));

        if( (abs(currentPtr->index(0) - endPtr->index(0)) <= termination_grid_num && 
             abs(currentPtr->index(1) - endPtr->index(1)) <= termination_grid_num && 
             abs(currentPtr->index(2) - endPtr->index(2)) <= termination_grid_num ) || shotHeu(currentPtr, endPtr) )
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
                        double tentative_fScore = current_gScore + expandPtr->edge_cost + getHeu(expandPtr, endPtr) * w_h;
                        double real_current_fScore = currentPtr->gScore + getHeu(currentPtr, endPtr) * w_h;

                        if(num_iter < 2)
                        {   
                            cout<<u(0)<<","<<u(1)<<","<<u(2)<<endl;
                            cout<<"current_gScore: "<<current_gScore<<endl;
                            cout<<"expandPtr->edge_cost: "<<expandPtr->edge_cost<<endl;
                            cout<<"Heu of expandPtr :"<<getHeu(expandPtr, endPtr)* w_h<<endl;
                            cout<<"tentative_fScore: "<<tentative_fScore<<endl;
                            cout<<"real_current_fScore fScore: "<<real_current_fScore + tie_breaker<<endl;
                        }

                        if(tentative_fScore <= real_current_fScore + tie_breaker)
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
                            neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
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
    if(is_shot_succ) // add shoting heuristic trajectory to the kino traj list
    {
        Vector3d coord;
        VectorXd poly1d, t;

        for(double time = resolution; time <= t_shot; time += resolution )
        {
            for ( int dim = 0; dim < 3; dim++ )
            {
                poly1d = coef_shot.row(dim);
                t = VectorXd::Zero(4);    

                for(int j = 0; j < 4; j ++)
                    t(j) = pow(time, j);

                coord(dim) = poly1d.dot(t);
            }

            statesList.push_back(coord);   
        }
    }

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

    succPtr->index = succIdx;
    succPtr->state = xt;
    succPtr->edge_cost = (u(0) * u(0) + u(1) * u(1) + u(2) * u(2) + w_t ) * t_exp;   

    auto finish = std::chrono::high_resolution_clock::now();
    time_in_forward += std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();                    

    return true;
}