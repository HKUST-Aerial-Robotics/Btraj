#include "fm_planer/trajectory_generator_lite.h"
using namespace std;    
using namespace Eigen;

#define inf 1>>30

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */


TrajectoryGenerator::TrajectoryGenerator(){}

TrajectoryGenerator::~TrajectoryGenerator(){}

MatrixXd TrajectoryGenerator::BezierPloyCoeffGeneration(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
            const VectorXd &C,
            const VectorXd &Cv,
            const VectorXd &Ca,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const int minimize_order,
            double & obj,
            const double margin)  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   

#if 1
#define MINORDER  minimize_order

    //ROS_WARN("[Bezier Trajectory] In TrajectoryGenerator");
    //ros::Time time_start = ros::Time::now();
    MatrixXd PolyCoeff;

    double initScale = corridor.front().t;
    double lstScale  = corridor.back().t;

    /*assert(_Time.size() == _Radius.size() );
    _Path = Path.block(1, 0, Path.rows() - 2, 3 );*/

    int _segment_num   = corridor.size();

    Vector3d StartPt  = pos.row(0); 
    Vector3d EndPt    = pos.row(1); 

    Vector3d StartVel = vel.row(0);
    Vector3d EndVel   = vel.row(1);

    Vector3d StartAcc = acc.row(0);
    Vector3d EndAcc   = acc.row(1);
    
//#####################################################################################################################
//Prepare for constaints and objective data stucture for Mosek solver .         
    // cout<<"maxVel: "<<maxVel<<endl;

    int n_poly = traj_order + 1;
    //cout<<"n_poly: "<<n_poly<<endl;
    int _s1d1CtrlP_num = n_poly;
  //  cout<<"_s1d1CtrlP_num: "<<_s1d1CtrlP_num<<endl;
    int _s1CtrlP_num   = 3 * _s1d1CtrlP_num;
/*    cout<<"_s1CtrlP_num: "<<_s1CtrlP_num<<endl;
    cout<<"_segment_num: "<<_segment_num<<endl;*/

    int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point
    int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point
    int equ_con_continuity_num = 3 * 3 * (_segment_num - 1);
    int _equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position
    
    int _inequ_con_num      = 0;
    int _high_order_con_num = 0; //3 * traj_order * _segment_num;
    int _CtrlP_num          = _segment_num * _s1CtrlP_num;
    

/*    cout<<"all variables number is: "<<_CtrlP_num<<endl;
    cout<<"all inequality number is: "<<_inequ_con_num<<endl;*/

    int con_num = _equ_con_num + _inequ_con_num + _high_order_con_num;
    //cout<<"all constraints number is: "<<con_num<<endl;

    double x_var[ _CtrlP_num];
    double primalobj;

    MSKrescodee  r; 
    /* ## define a container for constaints boundary and boundkey ## */ 
    /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > con_bdk; 
    
    //ROS_WARN("[Bezier Trajectory] Start stacking the bounding value");
    /***  Stack the bounding value for the linear inequality for the high order dynamic constraints  ***/
    // for velocity constraints

/*    for(int i = 0; i < _high_order_con_num; i++)
    {
        pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxVel,  + maxVel) );
        con_bdk.push_back(cb_ie);   
    }*/

    // for acceleration constraints
/*    for(int k =0; k < _segment_num; k++ ){
        for(int p = 0; p < 3; p++) // loop for x, y, z
            for(int i = 0; i < _s1d1CtrlP_num - 2; i++){
                pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - maxAcc,  maxAcc) ); 
                con_bdk.push_back(cb_ie);   
            }
    }*/

    //ROS_WARN("[Bezier Trajectory] equality bound %d", _equ_con_num);
    for(int i = 0; i < _equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = StartPt(i); 
        else if (i >= 3  && i < 6  ) beq_i = StartVel(i - 3); 
        else if (i >= 6  && i < 9  ) beq_i = StartAcc(i - 6);
        else if (i >= 9  && i < 12 ) beq_i = EndPt (i - 9 );
        else if (i >= 12 && i < 15 ) beq_i = EndVel(i - 12);
        else if (i >= 15 && i < 18 ) beq_i = EndAcc(i - 15);
        else beq_i = 0.0;

        //cout<<"beq_i: "<<beq_i<<endl;
        pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
        con_bdk.push_back(cb_eq);
    }

    /* ## define a container for control points' boundary and boundkey ## */ 
    /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
    vector< pair<MSKboundkeye, pair<double, double> > > var_bdk; 
    /*for(int i = 0; i < _CtrlP_num; i ++ )
    {
        pair<MSKboundkeye, pair<double, double> > vb_x  = make_pair( MSK_BK_FR, make_pair( - MSK_INFINITY, + MSK_INFINITY ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)
        var_bdk.push_back(vb_x);
    }*/ 
/*    double numerical_margin = 0.01;
    double eps = 0.00001;*/
    for(int k = 0; k < _segment_num; k++)
    {   
        Cube cube_     = corridor[k];
        double scale_k = cube_.t;

        for(int i = 0; i < 3; i++ )
        {   
            for(int j = 0; j < n_poly; j ++ )
            {   
                pair<MSKboundkeye, pair<double, double> > vb_x;

                double lo_bound = (cube_.box[i].first  + margin) / scale_k;
                double up_bound = (cube_.box[i].second - margin) / scale_k;

/*                if(fabs(up_bound) < eps)
                    up_bound = 0.0;

                if(fabs(lo_bound) < eps)
                    lo_bound = 0.0;*/
                //vb_x  = make_pair( MSK_BK_FR, make_pair(  - MSK_INFINITY, + MSK_INFINITY ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)
                vb_x  = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)

                var_bdk.push_back(vb_x);
            }
        } 
    }

    MSKint32t  j,i; 
    MSKenv_t   env; 
    MSKtask_t  task; 
    // Create the mosek environment. 
    r = MSK_makeenv( &env, NULL ); 
  
    // Create the optimization task. 
    r = MSK_maketask(env,con_num, _CtrlP_num, &task); 

// Parameters used in the optimizer
//######################################################################
    //MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_INTPNT );
    MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);
    MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-2);
    /*MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-3);
    MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-3 );*/
    
    //MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );

//######################################################################
    
    r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
    // Append 'con_num' empty constraints. 
     //The constraints will initially have no bounds. 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendcons(task,con_num);  

    // Append '_CtrlP_num' variables. The variables will initially be fixed at zero (x=0). 
    if ( r == MSK_RES_OK ) 
      r = MSK_appendvars(task,_CtrlP_num); 

    //ROS_WARN("set variables boundary");
    for(j = 0; j<_CtrlP_num && r == MSK_RES_OK; ++j){ 
      // Set the bounds on variable j : //  blx[j] <= x_j <= bux[j] 
        //cout<<"lower bound: "<<var_bdk[j].second.first<<" , upper bound: "<<var_bdk[j].second.second<<endl;
        if (r == MSK_RES_OK) 
            r = MSK_putvarbound(task, 
                                j,                            // Index of variable. 
                                var_bdk[j].first,             // Bound key.
                                var_bdk[j].second.first,      // Numerical value of lower bound.
                                var_bdk[j].second.second );   // Numerical value of upper bound.      
    } 
    
    //ROS_WARN("con_num: %d", con_num);
    //ROS_WARN("con_bdk size: %d", int(con_bdk.size() ));

    //ROS_WARN("set constraints boundary");
    // Set the bounds on constraints. 
    //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
    for( i = 0; i < con_num && r == MSK_RES_OK; i++ ) {
        //cout<<"bounding key: "<<con_bdk[i].first<<endl;
        //cout<<"lower bound: "<<con_bdk[i].second.first<<" , upper bound: "<<con_bdk[i].second.second<<endl;
        r = MSK_putconbound(task, 
                            i,                            // Index of constraint. 
                            con_bdk[i].first,             // Bound key.
                            con_bdk[i].second.first,      // Numerical value of lower bound.
                            con_bdk[i].second.second );   // Numerical value of upper bound. 
    }

    //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, inequality part");
    // Put variables in A by row one by one
    // #0   put the in-equality coinstraints in each control points
    // #0.1 positon constraints of each control point in the corresponding sphere
    // #0.2 velocity constraints of each 1st order control points within the limit
    // #0.3 acceleration constraints of each 2nd order control points within the limit
    int row_idx = 0;
    // The velocity constraints
/*    {   
        int sub_shift = 0;
        for(int k = 0; k < _segment_num ; k ++ )
        {   
            //int order = traj_order;
            for(int i = 0; i < 3; i++)
            {  // for x, y, z loop
                for(int p = 0; p < traj_order; p++)
                {
                    int nzi = 2;
                    MSKint32t asub[nzi];
                    double aval[nzi];

                    aval[0] = -1.0 * traj_order;
                    aval[1] =  1.0 * traj_order;

                    asub[0] = sub_shift + i * _s1d1CtrlP_num + p;    
                    asub[1] = sub_shift + i * _s1d1CtrlP_num + p + 1;    

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }
            }
            sub_shift += _s1CtrlP_num;
        }
    }*/

    // The acceleration constraints
/*    row_idx = row_idx_shift;
    for(int k = 0; k < _segment_num ; k ++ ){
        for(int i = 0; i < 3; i++){  // for x, y, z loop
            for(int p = 0; p < _s1d1CtrlP_num - 2; p++){
                
                int nzi = 3;
                MSKint32t asub[nzi];
                double aval[nzi];

                aval[0] =  1.0 * POLYORDER * (POLYORDER - 1) / _Scale(k);
                aval[1] = -2.0 * POLYORDER * (POLYORDER - 1) / _Scale(k);
                aval[2] =  1.0 * POLYORDER * (POLYORDER - 1) / _Scale(k);
                asub[0] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p;    
                asub[1] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p + 1;    
                asub[2] = k * _s1CtrlP_num + i * _s1d1CtrlP_num + p + 2;    
                
                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
        }
    }
    row_idx_shift += 3 * (_s1d1CtrlP_num - 2) * _segment_num;*/

     //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, equality part");   
    // #1   put the equality constraints in the start position
    // #1.1 positon constraints in the start point
    // #1.2 velocity constraints in the start point
    // #1.3 acceleration constraints in the start point
    // For position, velocity and acceleration, seperately
    
    ROS_WARN(" start position");
    /*   Start position  */
    {
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = 1.0 * initScale;
            asub[0] = i * _s1d1CtrlP_num;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = - 1.0 * traj_order;
            aval[1] =   1.0 * traj_order;
            asub[0] = i * _s1d1CtrlP_num;
            asub[1] = i * _s1d1CtrlP_num + 1;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);   
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] =   1.0 * traj_order * (traj_order - 1) / initScale;
            aval[1] = - 2.0 * traj_order * (traj_order - 1) / initScale;
            aval[2] =   1.0 * traj_order * (traj_order - 1) / initScale;
            asub[0] = i * _s1d1CtrlP_num;
            asub[1] = i * _s1d1CtrlP_num + 1;
            asub[2] = i * _s1d1CtrlP_num + 2;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }      
    // #2   put the equality constraints in the end position
    // #2.1 positon constraints in the end point
    // #2.2 velocity constraints in the end point
    // #2.3 acceleration constraints in the end point

    /*   End position  */
    //ROS_WARN(" end position");
    {   
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            //cout<<"asub[0]: "<<_CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num<<endl;
            asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
            aval[0] = 1.0 * lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 1;
            asub[1] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
            aval[0] = - 1.0;
            aval[1] =   1.0;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 2;
            asub[1] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num - 1;
            asub[2] = _CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num;
            aval[0] =   1.0 / lstScale;
            aval[1] = - 2.0 / lstScale;
            aval[2] =   1.0 / lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
    }

    // #3   put the equality coinstraints in each joint positions 
    // #3.1 positon constraints in each joint positions
    // #3.2 velocity constraints in each joint positions
    // #3.3 acceleration constraints in each joint positions

    /*   joint points  */
    //ROS_WARN(" joint position");
    {
        int sub_shift = 0;
        double val0, val1;
        for(int k = 0; k < (_segment_num - 1); k ++ )
        {   
            double scale_k = corridor[k].t;
            double scale_n = corridor[k+1].t;
            // position :
            val0 = scale_k;
            val1 = scale_n;
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 2;
                MSKint32t asub[nzi];
                double aval[nzi];

                // This segment's last control point
                aval[0] = 1.0 * val0;
                asub[0] = sub_shift + (i+1) * _s1d1CtrlP_num - 1;
                /*cout<<"idx 0: "<<sub_shift + (i+1) * _s1d1CtrlP_num - 1<<endl;
                cout<<"'row_idx: "<<row_idx<<endl;*/
                //asub[0] = k * _s1CtrlP_num + (i+1) * _s1d1CtrlP_num - 1;
                // Next segment's first control point
                aval[1] = -1.0 * val1;
                asub[1] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num;
                //cout<<"idx 1: "<<sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num<<endl;
                //asub[1] = (k+1) * _s1CtrlP_num + i * _s1d1CtrlP_num;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            // velocity :
          /*  val0 = 1.0;
            val1 = 1.0;*/
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 4;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] = -1.0;
                aval[1] =  1.0;
               /* cout<<"idx 0: "<<sub_shift + (i+1) * _s1d1CtrlP_num - 2<<endl;
                cout<<"idx 1: "<<sub_shift + (i+1) * _s1d1CtrlP_num - 1<<endl;*/
                asub[0] = sub_shift + (i+1) * _s1d1CtrlP_num - 2;    
                asub[1] = sub_shift + (i+1) * _s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[2] =  1.0;
                aval[3] = -1.0;

                /*cout<<"idx 2: "<<sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num<<endl;
                cout<<"idx 3: "<<sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num + 1<<endl;*/

                asub[2] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num;    
                asub[3] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num + 1;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }
            // acceleration :
            val0 = 1.0 / scale_k;
            val1 = 1.0 / scale_n;
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 6;
                MSKint32t asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] =  1.0  * val0;
                aval[1] = -2.0  * val0;
                aval[2] =  1.0  * val0;
                asub[0] = sub_shift + (i+1) * _s1d1CtrlP_num - 3;    
                asub[1] = sub_shift + (i+1) * _s1d1CtrlP_num - 2;   
                asub[2] = sub_shift + (i+1) * _s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[3] =  -1.0  * val1;
                aval[4] =   2.0  * val1;
                aval[5] =  -1.0  * val1;
                asub[3] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num;    
                asub[4] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num + 1;
                asub[5] = sub_shift + _s1CtrlP_num + i * _s1d1CtrlP_num + 2;

                r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                row_idx ++;
            }

            sub_shift += _s1CtrlP_num;
        }
    }

    //ROS_WARN("[Bezier Trajectory] Start stacking the objective");
    //cout<<"MQM: \n"<<MQM<<endl;
    // ROS_WARN("[Solver] prepare for stacking the objective");
    
    int NUMQNZ = 0;
    for(int i = 0; i < _segment_num; i ++)
    {   
        int NUMQ_blk = (traj_order + 1);                       // default minimize the jerk and MINORDER = 3
        NUMQNZ      += 3 * NUMQ_blk * (NUMQ_blk + 1) / 2;
    }
    MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
    double     qval[NUMQNZ];
    
    {    
        int sub_shift = 0;
        int idx = 0;
        for(int k = 0; k < _segment_num; k ++)
        {
            double scale_k = corridor[k].t;
            for(int p = 0; p < 3; p ++ )
                for( int i = 0; i < _s1d1CtrlP_num; i ++ )
                    for( int j = 0; j < _s1d1CtrlP_num; j ++ )
                        if( i >= j )
                        {
                            qsubi[idx] = sub_shift + p * _s1d1CtrlP_num + i;   
                            qsubj[idx] = sub_shift + p * _s1d1CtrlP_num + j;  
                            qval[idx]  = MQM(i, j) /(double)pow(scale_k, 3);// + 0.1;// + 0.1;                
                            idx ++ ;
                        }

            sub_shift += _s1CtrlP_num;
        }
    } 
    
    ros::Time time_end1 = ros::Time::now();
    //ROS_WARN("Time in variables stack is");
    //cout<<time_end1 - time_start<<endl;

    if ( r== MSK_RES_OK )
         r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); 
    
    if ( r==MSK_RES_OK ) 
         r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    
    //ros::Time time_opt = ros::Time::now();
    bool solve_ok = false;
    if ( r==MSK_RES_OK ) 
      { 
        //ROS_WARN("Prepare to solve the problem ");   
        MSKrescodee trmcode; 
        r = MSK_optimizetrm(task,&trmcode); 
        MSK_solutionsummary (task,MSK_STREAM_LOG); 
          
        if ( r==MSK_RES_OK ) 
        { 
          MSKsolstae solsta; 
          MSK_getsolsta (task,MSK_SOL_ITR,&solsta); 
           
          switch(solsta) 
          { 
            case MSK_SOL_STA_OPTIMAL:    
            case MSK_SOL_STA_NEAR_OPTIMAL: 
              
            
            r = MSK_getxx(task, 
                          MSK_SOL_ITR,    // Request the interior solution.  
                          x_var); 
            
            r = MSK_getprimalobj(
                task,
                MSK_SOL_ITR,
                &primalobj);

            obj = primalobj;
            solve_ok = true;
            
            break; 
            
            case MSK_SOL_STA_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_PRIM_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:   
              printf("Primal or dual infeasibility certificate found.\n"); 
              break; 
               
            case MSK_SOL_STA_UNKNOWN: 
              printf("The status of the solution could not be determined.\n"); 
              //solve_ok = true; // debug
              break; 
            default: 
              printf("Other solution status."); 
              break; 
          } 
        } 
        else 
        { 
          printf("Error while optimizing.\n"); 
        } 
      }
     
      if (r != MSK_RES_OK) 
      { 
        // In case of an error print error code and description. 
        char symname[MSK_MAX_STR_LEN]; 
        char desc[MSK_MAX_STR_LEN]; 
         
        printf("An error occurred while optimizing.\n");      
        MSK_getcodedesc (r, 
                         symname, 
                         desc); 
        printf("Error %s - '%s'\n",symname,desc); 
      } 
    
    MSK_deletetask(&task); 
    MSK_deleteenv(&env); 

    ros::Time time_end2 = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_end2 - time_end1<<endl;

    if(!solve_ok){
      MatrixXd poly_fail = MatrixXd::Identity(3,3);
      ROS_WARN("In solver, falied ");
      return poly_fail;
    }

    VectorXd d_var(_CtrlP_num);
    for(int i = 0; i < _CtrlP_num; i++)
        d_var(i) = x_var[i];

    //cout<<"solution: "<<d_var<<endl;

    //ROS_BREAK();
    
    PolyCoeff = MatrixXd::Zero(_segment_num, 3 *(traj_order + 1) );

    int var_shift = 0;
    for(int i = 0; i < _segment_num; i++ )
    {
        for(int j = 0; j < 3 * n_poly; j++)
            PolyCoeff(i , j) = d_var(j + var_shift);

        var_shift += 3 * n_poly;
    }   

    //cout<<"PolyCoeff:\n"<<PolyCoeff<<endl;

    return PolyCoeff;
#endif
}

MatrixXd TrajectoryGenerator::BezierPloyCoeffGenerationLite(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
            const VectorXd &C,
            const VectorXd &Cv,
            const VectorXd &Ca,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const int minimize_order,
            double & obj,
            const double margin ) // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   

#if 1
#define MINORDER  minimize_order

    //ROS_WARN("[Bezier Trajectory] In TrajectoryGenerator");
    //ros::Time time_start = ros::Time::now();
    MatrixXd PolyCoeff;

    double initScale = corridor.front().t;
    double lstScale  = corridor.back().t;

    cout<<"initScale: "<<initScale<<endl;
    cout<<"lstScale:  "<<lstScale<<endl;
    /*assert(_Time.size() == _Radius.size() );

    _Path = Path.block(1, 0, Path.rows() - 2, 3 );*/
        
    int _segment_num   = corridor.size();

    cout<<"_segment_num: "<<_segment_num<<endl;

    Vector3d StartPt  = pos.row(0); 
    Vector3d EndPt    = pos.row(1); 

    Vector3d StartVel = vel.row(0);
    Vector3d EndVel   = vel.row(1);

    Vector3d StartAcc = acc.row(0);
    Vector3d EndAcc   = acc.row(1);

    int n_poly = traj_order + 1;
    int _s1d1CtrlP_num = n_poly;
    int _s1CtrlP_num   = _s1d1CtrlP_num;

    int equ_con_s_num = 3; // p, v, a in x, y, z axis at the start point
    int equ_con_e_num = 3; // p, v, a in x, y, z axis at the end point
    int equ_con_continuity_num = 3 * (_segment_num - 1);
    int _equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position
    
    int _inequ_con_num      = 0;
    int _high_order_con_num = traj_order * _segment_num;
    int _CtrlP_num          = _segment_num * _s1CtrlP_num;

    int con_num = _equ_con_num + _inequ_con_num + _high_order_con_num;
    double x_var[ _CtrlP_num];
    
    double primalobj;
    MSKrescodee  r; 
    /* ## define a container for constaints boundary and boundkey ## */ 
    /* ## dataType in the double pair is : boundary type, lower bound, upper bound ## */

    ros::Time time_end1 = ros::Time::now();
    MatrixXd d_var(3, _CtrlP_num);
    obj = 0.0;
    for(int p = 0; p < 3; p++)
    {
        vector< pair<MSKboundkeye, pair<double, double> > > con_bdk;
        
        for(int i = 0; i < _high_order_con_num; i++)
        {
            pair<MSKboundkeye, pair<double, double> > cb_ie = make_pair( MSK_BK_RA, make_pair( - 1000.0,  + 1000.0) );
            con_bdk.push_back(cb_ie);   
        }

        ROS_WARN("[Bezier Trajectory] equality bound %d", _equ_con_num);
        for(int i = 0; i < _equ_con_num; i ++ )
        { 
            double beq_i;
            if(i < 1)                    beq_i = StartPt(p); 
            else if (i >= 1  && i < 2  ) beq_i = StartVel(p); 
            else if (i >= 2  && i < 3  ) beq_i = StartAcc(p);
            else if (i >= 3  && i < 4 ) beq_i = EndPt (p);
            else if (i >= 4  && i < 5 ) beq_i = EndVel(p);
            else if (i >= 5  && i < 6 ) beq_i = EndAcc(p);
            else beq_i = 0.0;

            //cout<<"beq_i: "<<beq_i<<endl;
            pair<MSKboundkeye, pair<double, double> > cb_eq = make_pair( MSK_BK_FX, make_pair( beq_i, beq_i ) ); // # cb_eq means: constriants boundary of equality constrain
            con_bdk.push_back(cb_eq);
        }

        /* ## define a container for control points' boundary and boundkey ## */ 
        /* ## dataType in one tuple is : boundary type, lower bound, upper bound ## */
        vector< pair<MSKboundkeye, pair<double, double> > > var_bdk;
        double eps = 0.00001;
        for(int k = 0; k < _segment_num; k++)
        {   
            Cube cube_     = corridor[k];
            double scale_k = cube_.t;

            for(int j = 0; j < n_poly; j ++ )
            {   
                pair<MSKboundkeye, pair<double, double> > vb_x;

                double lo_bound = (cube_.box[p].first  + margin) / scale_k;
                double up_bound = (cube_.box[p].second - margin) / scale_k;

                if(fabs(up_bound) < eps)
                    up_bound = 0.0;

                if(fabs(lo_bound) < eps)
                    lo_bound = 0.0;
                //vb_x  = make_pair( MSK_BK_FR, make_pair(  - MSK_INFINITY, + MSK_INFINITY ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)
                vb_x  = make_pair( MSK_BK_RA, make_pair( lo_bound, up_bound ) ); // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)

                var_bdk.push_back(vb_x);
            }
        }

        MSKint32t  j,i; 
        MSKenv_t   env; 
        MSKtask_t  task; 
        // Create the mosek environment. 
        r = MSK_makeenv( &env, NULL ); 
      
        // Create the optimization task. 
        r = MSK_maketask(env,con_num, _CtrlP_num, &task); 

    // Parameters used in the optimizer
    //######################################################################
        MSK_putintparam (task, MSK_IPAR_OPTIMIZER , MSK_OPTIMIZER_INTPNT );
        MSK_putintparam (task, MSK_IPAR_NUM_THREADS, 1);
        MSK_putdouparam (task, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 1e-2);
        
        MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-3);
        MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-3);
        MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-3 );
        
        //MSK_putdouparam (task, MSK_DPAR_INTPNT_TOL_REL_GAP, 5e-2 );

    //######################################################################
        
        r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr); 
        // Append 'con_num' empty constraints. 
         //The constraints will initially have no bounds. 
        if ( r == MSK_RES_OK ) 
          r = MSK_appendcons(task,con_num);  

        // Append '_CtrlP_num' variables. The variables will initially be fixed at zero (x=0). 
        if ( r == MSK_RES_OK ) 
          r = MSK_appendvars(task,_CtrlP_num); 

        ROS_WARN("set variables boundary");
        for(j = 0; j<_CtrlP_num && r == MSK_RES_OK; ++j)
        { 
            // Set the bounds on variable j : //  blx[j] <= x_j <= bux[j] 
            // cout<<"lower bound: "<<var_bdk[j].second.first<<" , upper bound: "<<var_bdk[j].second.second<<endl;
            if (r == MSK_RES_OK) 
                r = MSK_putvarbound(task, 
                                    j,                            // Index of variable. 
                                    var_bdk[j].first,             // Bound key.
                                    var_bdk[j].second.first,      // Numerical value of lower bound.
                                    var_bdk[j].second.second );   // Numerical value of upper bound.      
        } 
        
        ROS_WARN("con_num: %d", con_num);
        ROS_WARN("con_bdk size: %d", int(con_bdk.size() ));

        ROS_WARN("set constraints boundary");
        // Set the bounds on constraints. 
        //   for i=1, ...,con_num : blc[i] <= constraint i <= buc[i] 
        for( i = 0; i < con_num && r == MSK_RES_OK; i++ ) {
            //cout<<"bounding key: "<<con_bdk[i].first<<endl;
            //cout<<"lower bound: "<<con_bdk[i].second.first<<" , upper bound: "<<con_bdk[i].second.second<<endl;
            r = MSK_putconbound(task, 
                                i,                            // Index of constraint. 
                                con_bdk[i].first,             // Bound key.
                                con_bdk[i].second.first,      // Numerical value of lower bound.
                                con_bdk[i].second.second );   // Numerical value of upper bound. 
        }

        //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, inequality part");
        // Put variables in A by row one by one
        // #0   put the in-equality coinstraints in each control points
        // #0.1 positon constraints of each control point in the corresponding sphere
        // #0.2 velocity constraints of each 1st order control points within the limit
        // #0.3 acceleration constraints of each 2nd order control points within the limit
        int row_idx = 0;
        
        // The velocity constraints
        {   
            int sub_shift = 0;
            for(int k = 0; k < _segment_num ; k ++ )
            {   
                //int order = traj_order;
                for(int p = 0; p < traj_order; p++)
                {
                    int nzi = 2;
                    MSKint32t asub[nzi];
                    double aval[nzi];

                    aval[0] = -1.0 * traj_order;
                    aval[1] =  1.0 * traj_order;

                    asub[0] = sub_shift + p;    
                    asub[1] = sub_shift + p + 1;    

                    /*cout<<"val: "<<endl;
                    cout<<aval[0]<<endl;
                    cout<<aval[1]<<endl;
                    cout<<"sub: "<<endl;
                    cout<<asub[0]<<endl;
                    cout<<asub[1]<<endl;*/

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }

                sub_shift += _s1CtrlP_num;
            }
        }

        ROS_WARN(" start position");
        /*   Start position  */
        {
            // position :
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = 1.0 * initScale;
            asub[0] = 0;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
            // velocity :
        {
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] = - 1.0 * traj_order;
            aval[1] =   1.0 * traj_order;
            asub[0] = 0;
            asub[1] = 1;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);   
            row_idx ++;
        }            
            // acceleration : 
        {
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            aval[0] =   1.0 * traj_order * (traj_order - 1) / initScale;
            aval[1] = - 2.0 * traj_order * (traj_order - 1) / initScale;
            aval[2] =   1.0 * traj_order * (traj_order - 1) / initScale;
            asub[0] = 0;
            asub[1] = 1;
            asub[2] = 2;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }      

        // #2   put the equality constraints in the end position
        // #2.1 positon constraints in the end point
        // #2.2 velocity constraints in the end point
        // #2.3 acceleration constraints in the end point
        /*   End position  */
        ROS_WARN(" end position");
        {   
            // position :
            int nzi = 1;
            MSKint32t asub[nzi];
            double aval[nzi];
            //cout<<"asub[0]: "<<_CtrlP_num - 1 - (2 - i) * _s1d1CtrlP_num<<endl;
            asub[0] = _CtrlP_num - 1;
            aval[0] = 1.0 * lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }
            // velocity :
        {
            int nzi = 2;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = _CtrlP_num - 1 - 1;
            asub[1] = _CtrlP_num - 1;
            aval[0] = - 1.0;
            aval[1] =   1.0;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }    
            // acceleration : 
        {
            int nzi = 3;
            MSKint32t asub[nzi];
            double aval[nzi];
            asub[0] = _CtrlP_num - 1 - 2;
            asub[1] = _CtrlP_num - 1 - 1;
            asub[2] = _CtrlP_num - 1;
            aval[0] =   1.0 / lstScale;
            aval[1] = - 2.0 / lstScale;
            aval[2] =   1.0 / lstScale;
            r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            row_idx ++;
        }

        // #3   put the equality coinstraints in each joint positions 
        // #3.1 positon constraints in each joint positions
        // #3.2 velocity constraints in each joint positions
        // #3.3 acceleration constraints in each joint positions

        /*   joint points  */
        ROS_WARN(" joint position");
        {
            int sub_shift = 0;
            double val0, val1;
            for(int k = 0; k < (_segment_num - 1); k ++ )
            {   
                double scale_k = corridor[k].t;
                double scale_n = corridor[k+1].t;
                // position :
                {
                    val0 = scale_k;
                    val1 = scale_n;

                    int nzi = 2;
                    MSKint32t asub[nzi];
                    double aval[nzi];

                    // This segment's last control point
                    aval[0] = 1.0 * val0;
                    asub[0] = sub_shift + _s1d1CtrlP_num - 1;

                    aval[1] = -1.0 * val1;
                    asub[1] = sub_shift + _s1CtrlP_num;
                    
                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }
                // velocity :
                {
                    int nzi = 4;
                    MSKint32t asub[nzi];
                    double aval[nzi];
                    
                    // This segment's last velocity control point
                    aval[0] = -1.0;
                    aval[1] =  1.0;

                    asub[0] = sub_shift + _s1d1CtrlP_num - 2;    
                    asub[1] = sub_shift + _s1d1CtrlP_num - 1;   
                    // Next segment's first velocity control point
                    aval[2] =  1.0;
                    aval[3] = -1.0;

                    asub[2] = sub_shift + _s1CtrlP_num;    
                    asub[3] = sub_shift + _s1CtrlP_num + 1;

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }
                // acceleration :
                {
                    val0 = 1.0 / scale_k;
                    val1 = 1.0 / scale_n;

                    int nzi = 6;
                    MSKint32t asub[nzi];
                    double aval[nzi];
                    
                    // This segment's last velocity control point
                    aval[0] =  1.0  * val0;
                    aval[1] = -2.0  * val0;
                    aval[2] =  1.0  * val0;
                    asub[0] = sub_shift + _s1d1CtrlP_num - 3;    
                    asub[1] = sub_shift + _s1d1CtrlP_num - 2;   
                    asub[2] = sub_shift + _s1d1CtrlP_num - 1;   
                    // Next segment's first velocity control point
                    aval[3] =  -1.0  * val1;
                    aval[4] =   2.0  * val1;
                    aval[5] =  -1.0  * val1;
                    asub[3] = sub_shift + _s1CtrlP_num;    
                    asub[4] = sub_shift + _s1CtrlP_num + 1;
                    asub[5] = sub_shift + _s1CtrlP_num + 2;

                    r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                    row_idx ++;
                }

                sub_shift += _s1CtrlP_num;
            }
        }

        ROS_WARN("[Bezier Trajectory] Start stacking the objective");
        //cout<<"MQM: \n"<<MQM<<endl;
        // ROS_WARN("[Solver] prepare for stacking the objective");
        
        int NUMQNZ = 0;
        for(int i = 0; i < _segment_num; i ++)
        {   
            int NUMQ_blk = (traj_order + 1);                       // default minimize the jerk and MINORDER = 3
            NUMQNZ      += NUMQ_blk * (NUMQ_blk + 1) / 2;
        }
        MSKint32t  qsubi[NUMQNZ], qsubj[NUMQNZ];
        double     qval[NUMQNZ];
        
        {    
            int sub_shift = 0;
            int idx = 0;
            for(int k = 0; k < _segment_num; k ++)
            {
                double scale_k = corridor[k].t;
                for( int i = 0; i < _s1d1CtrlP_num; i ++ )
                    for( int j = 0; j < _s1d1CtrlP_num; j ++ )
                        if( i >= j )
                        {
                            qsubi[idx] = sub_shift + i;   
                            qsubj[idx] = sub_shift + j;  
                            qval[idx]  = MQM(i, j) / pow(scale_k, 3);// + 0.1;// + 0.1;                
                            idx ++ ;
                        }

                sub_shift += _s1CtrlP_num;
            }
        } 
        
        ros::Time time_end1 = ros::Time::now();
        //ROS_WARN("Time in variables stack is");
        //cout<<time_end1 - time_start<<endl;

        if ( r== MSK_RES_OK )
             r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval); 
        
        if ( r==MSK_RES_OK ) 
             r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
        
        //ros::Time time_opt = ros::Time::now();
        bool solve_ok = false;
        if ( r==MSK_RES_OK ) 
          { 
            //ROS_WARN("Prepare to solve the problem ");   
            MSKrescodee trmcode; 
            r = MSK_optimizetrm(task,&trmcode); 
            MSK_solutionsummary (task,MSK_STREAM_LOG); 
              
            if ( r==MSK_RES_OK ) 
            { 
              MSKsolstae solsta; 
              MSK_getsolsta (task,MSK_SOL_ITR,&solsta); 
               
              switch(solsta) 
              { 
                case MSK_SOL_STA_OPTIMAL:    
                case MSK_SOL_STA_NEAR_OPTIMAL: 
                  
                
                r = MSK_getxx(task, 
                              MSK_SOL_ITR,    // Request the interior solution.  
                              x_var); 

                r = MSK_getprimalobj(
                task,
                MSK_SOL_ITR,
                &primalobj);

                obj += primalobj;

                solve_ok = true;
                
                break; 
                
                case MSK_SOL_STA_DUAL_INFEAS_CER: 
                case MSK_SOL_STA_PRIM_INFEAS_CER: 
                case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
                case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:   
                  printf("Primal or dual infeasibility certificate found.\n"); 
                  break; 
                   
                case MSK_SOL_STA_UNKNOWN: 
                  printf("The status of the solution could not be determined.\n"); 
                  //solve_ok = true; // debug
                  break; 
                default: 
                  printf("Other solution status."); 
                  break; 
              } 
            } 
            else 
            { 
              printf("Error while optimizing.\n"); 
            } 
          }
         
          if (r != MSK_RES_OK) 
          { 
            // In case of an error print error code and description. 
            char symname[MSK_MAX_STR_LEN]; 
            char desc[MSK_MAX_STR_LEN]; 
             
            printf("An error occurred while optimizing.\n");      
            MSK_getcodedesc (r, 
                             symname, 
                             desc); 
            printf("Error %s - '%s'\n",symname,desc); 
          } 
        
        MSK_deletetask(&task); 
        MSK_deleteenv(&env); 
        
        if(!solve_ok){
          MatrixXd poly_fail = MatrixXd::Identity(3,3);
          ROS_WARN("In solver, falied ");
          return poly_fail;
        }
        
        for(int i = 0; i < _CtrlP_num; i++)
            d_var(p, i) = x_var[i];

       /* ROS_WARN("cehck result in one dimension");
        cout<<d_var.row(p)<<endl;*/
    }

    ros::Time time_end2 = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_end2 - time_end1<<endl;

    PolyCoeff = MatrixXd::Zero(_segment_num, 3 *(traj_order + 1) );

    int var_shift = 0;
    for(int k = 0; k < _segment_num; k++ )
    {   
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < n_poly; j++)
                PolyCoeff(k , i * _s1CtrlP_num + j) = d_var(i, j + var_shift);
            
        var_shift += _s1CtrlP_num;
    }   

    //cout<<"PolyCoeff:\n"<<PolyCoeff<<endl;

    return PolyCoeff;
#endif
}