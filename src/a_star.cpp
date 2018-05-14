#include "pathFinding.h"

void gridPathFinder::initGridNodeMap(double _resolution, Vector3d global_xyz_l)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    cout<<"gl_xl: "<<gl_xl<<endl;
    cout<<"gl_yl: "<<gl_yl<<endl;
    cout<<"gl_zl: "<<gl_zl<<endl;
    
    resolution = _resolution;
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++){
       GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
       for(int j = 0; j < GLY_SIZE; j++){
           GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
           for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
           }
       }
    }
}

void gridPathFinder::linkLocalMap(int *** map, Vector3d xyz_l)
{
    gridMap = map;
    loc_map_o = xyz_l;
    Vector3i idx_o = coord2gridIndex(loc_map_o);
    for(int64_t i = 0; i < X_SIZE; i++){
        for(int64_t j = 0; j < Y_SIZE; j++){
            for(int64_t k = 0; k < Z_SIZE; k++){
                GridNodePtr ptr = GridNodeMap[i+idx_o(0)][j+idx_o(1)][k+idx_o(2)];
                ptr->id = 0;
                ptr->occupancy = map[i][j][k];
            }
        }
    }
}

void gridPathFinder::resetLocalMap()
{   
//    Vector3i idx_o = coord2gridIndex(loc_map_o);
    //ROS_WARN("grid map reset");
    //cout<<"local map origin: \n"<<idx_o<<endl;

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
        GridNodePtr tmpPtr = ptr.second;
        tmpPtr->occupancy = 0; // forget the occupancy
        tmpPtr->id = 0;
        tmpPtr->cameFrom = NULL;
        tmpPtr->gScore = inf;
        tmpPtr->fScore = inf;
    }

    expandedNodes.clear();

/*    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++)
            {   
                GridNodePtr tmpPtr = GridNodeMap[i][j][k];
                
                if(tmpPtr->id == -1)
                    cout<<i<<" , "<<j<<" , "<<k<<endl;

                tmpPtr->occupancy = 0; // forget the occupancy
                tmpPtr->id = 0;
                tmpPtr->cameFrom = NULL;
                tmpPtr->gScore = inf;
                tmpPtr->fScore = inf;
            }*/

    //ROS_WARN("local map reset finish");
}

GridNodePtr gridPathFinder::pos2gridNodePtr(Vector3d pos)
{
    Vector3i idx = coord2gridIndex(pos);
    GridNodePtr grid_ptr = new GridNode(idx, pos);

    return grid_ptr;
}

Vector3d gridPathFinder::gridIndex2coord(Vector3i index)
{
  Vector3d pt;
  pt(0) = index(0) * resolution + gl_xl + 0.5 * resolution;
  pt(1) = index(1) * resolution + gl_yl + 0.5 * resolution;
  pt(2) = index(2) * resolution + gl_zl + 0.5 * resolution;
  return pt;
}

Vector3i gridPathFinder::coord2gridIndex(Vector3d pt)
{
  Vector3i idx;
  idx << min(max( int((pt(0) - gl_xl) / resolution), 0), GLX_SIZE - 1),
         min(max( int((pt(1) - gl_yl) / resolution), 0), GLY_SIZE - 1),
         min(max( int((pt(2) - gl_zl) / resolution), 0), GLZ_SIZE - 1);      

  return idx;
}

double gridPathFinder::getDiagHeu(GridNodePtr node1, GridNodePtr node2)
{   
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    double h;
    int diag = min(min(dx, dy), dz);
    dx -= diag;
    dy -= diag;
    dz -= diag;

    if (dx == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
    }
    if (dy == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
    }
    if (dz == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
    }
    return h;
}

double gridPathFinder::getManhHeu(GridNodePtr node1, GridNodePtr node2)
{   
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

double gridPathFinder::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{   
    return (node2->index - node1->index).norm();
}

double gridPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    return tie_breaker * getDiagHeu(node1, node2);
    //return tie_breaker * getEuclHeu(node1, node2);
}

vector<GridNodePtr> gridPathFinder::retrievePath(GridNodePtr current)
{   
    vector<GridNodePtr> path;
    path.push_back(current);

    while(current->cameFrom != NULL)
    {
        current = current -> cameFrom;
        path.push_back(current);
    }

    return path;
}

vector<GridNodePtr> gridPathFinder::getVisitedNodes()
{   
    vector<GridNodePtr> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++)
            {   
                if(GridNodeMap[i][j][k]->id != 0)
                //if(GridNodeMap[i][j][k]->id == -1)
                    visited_nodes.push_back(GridNodeMap[i][j][k]);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

void gridPathFinder::AstarSearch(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt)
{   
    ros::Time time_1 = ros::Time::now();    
    //updateGridMap();
    /*ros::Time time_m = ros::Time::now();
    ROS_WARN(" Time consume in stack variables for A star is %f", (time_m - time_1).toSec());
    ROS_WARN("[A star] GridNodeMap initialized...");
*/
/*    cout<<"gl_xl: "<<gl_xl<<endl;
    cout<<"gl_yl: "<<gl_yl<<endl;
    cout<<"gl_zl: "<<gl_zl<<endl;

    cout<<"GLX_SIZE: "<<GLX_SIZE<<endl;
    cout<<"GLY_SIZE: "<<GLY_SIZE<<endl;
    cout<<"GLZ_SIZE: "<<GLZ_SIZE<<endl;*/

    GridNodePtr startPtr = pos2gridNodePtr(start_pt);
    GridNodePtr endPtr   = pos2gridNodePtr(end_pt);

    openSet.clear();
/*
    ROS_WARN("[Astar]path finder starts to find path:");
    cout << "from: " << endl << start_pt << endl << "to: " << endl << end_pt << endl;

    ROS_WARN("[Astar]path finder actually find path after rounding: "); 
    cout << "from: " << endl << startPtr->real_coord << endl << "to: " << endl << endPtr->real_coord << endl;       */

    GridNodePtr neighborPtr = NULL;
    GridNodePtr current = NULL;

    startPtr -> gScore = 0;
    startPtr -> fScore = getHeu(startPtr, endPtr);
    startPtr -> id = 1; //put start node in open set
    startPtr -> coord = start_pt;
    openSet.insert( make_pair(startPtr -> fScore, startPtr) ); //put start in open set

    double tentative_gScore;

    int num_iter = 0;
    while ( !openSet.empty() )
    {   
        num_iter ++;
        current = openSet.begin() -> second;

        if(current->index(0) == endPtr->index(0)
        && current->index(1) == endPtr->index(1)
        && current->index(2) == endPtr->index(2) )
        {
            ROS_WARN("[Astar]Reach goal..");
            //cout << "goal coord: " << endl << current->real_coord << endl; 
            cout << "total number of iteration used in Astar: " << num_iter  << endl;
            ros::Time time_2 = ros::Time::now();
            ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
            gridPath = retrievePath(current);
            return;
        }         
        openSet.erase(openSet.begin());
        current -> id = -1; //move current node from open set to closed set.
        expandedNodes.push_back(current);

        for(int dx = -1; dx < 2; dx++)
            for(int dy = -1; dy < 2; dy++)
                for(int dz = -1; dz < 2; dz++){
                    if(dx == 0 && dy == 0 && dz ==0){
                        continue; //skip the current point itself
                    }
                    /*if(abs(dx) == 1 && abs(dy) == 1 && abs(dz) ==1){
                        continue; //skip the current point itself
                    }*/

                    Vector3i neighborIdx;
                    neighborIdx(0) = (current -> index)(0) + dx;
                    neighborIdx(1) = (current -> index)(1) + dy;
                    neighborIdx(2) = (current -> index)(2) + dz;

                    if(    neighborIdx(0) < 0 || neighborIdx(0) >= GLX_SIZE
                        || neighborIdx(1) < 0 || neighborIdx(1) >= GLY_SIZE
                        || neighborIdx(2) < 0 || neighborIdx(2) >= GLZ_SIZE){
                        continue;
                    }

                    neighborPtr = GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];

                    if(neighborPtr->occupancy > 0.5){
                        continue;
                    }

                    if(neighborPtr -> id == -1){
                        continue; //ignore neighbor which is already evaluated (in closed set).
                    }

                    double static_cost = sqrt(dx * dx + dy * dy + dz * dz);
                    
                    tentative_gScore = current -> gScore + static_cost; 

                    if(neighborPtr -> id != 1){
                        //discover a new node
                        neighborPtr -> id        = 1;
                        neighborPtr -> cameFrom  = current;
                        neighborPtr -> gScore    = tentative_gScore;
                        neighborPtr -> fScore    = neighborPtr -> gScore + getHeu(neighborPtr, endPtr); 
                        neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                        continue;
                    }
                    else if(tentative_gScore <= neighborPtr-> gScore){ //in open set and need update
                        neighborPtr -> cameFrom = current;
                        neighborPtr -> gScore = tentative_gScore;
                        neighborPtr -> fScore = tentative_gScore + getHeu(neighborPtr, endPtr); 
                        openSet.erase(neighborPtr -> nodeMapIt);
                        neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                    }
                        
                }
    }

    ros::Time time_2 = ros::Time::now();
    ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
}

vector<GridNodePtr> gridPathFinder::getPath()
{
    return gridPath;
}

void gridPathFinder::resetPath()
{
    gridPath.clear();
}

bool gridPathFinder::CheckGuidePathCollision()
{   
    if(gridPath.size() == 0) // no path exists
        return true;

    GridNodePtr nodePtr = NULL;
    for(int i = 0; i < (int)gridPath.size(); i++)
    {
        nodePtr = gridPath[i];
        if(nodePtr->occupancy == 1)
            return true;
    }

    return false;
}