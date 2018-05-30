#include "bezier_base.h"

MatrixXd Bernstein::CholeskyDecomp(MatrixXd Q) // return square root F of Q; Q = F' * F
{
	MatrixXd F, Ft;
	Eigen::LDLT< MatrixXd > ldlt(Q);
    F = ldlt.matrixL();
    F = ldlt.transpositionsP().transpose() * F;
    F *= ldlt.vectorD().array().sqrt().matrix().asDiagonal();
	Ft = F.transpose();

	return Ft;
}

int Bernstein::setParam(int poly_order_min, int poly_order_max, double min_order)
{
	int ret = (poly_order_min >=3 && poly_order_max <= 13) ? 1 : -1; // Now we only applied for poly order from 3 ~ 12 ( num of control from 4 ~ 13 )

	_order_min = poly_order_min;
	_order_max = poly_order_max;
	_min_order = min_order; 

	const static auto factorial = [](int n)
	{
        int fact = 1;

        for(int i = n; i > 0 ; i--)
          fact *= i;

        return fact; 
	};

	const static auto combinatorial = [](int n, int k) // for calculate n choose k combination problem
	{
		return factorial(n) / (factorial(k) * factorial(n - k));
	};

	CList.clear();
	CvList.clear();
	CaList.clear();
	CjList.clear();
	
	for(int order = 0; order <= _order_max; order++)
	{	
		MatrixXd M;   // Mapping matrix, used to map the coefficients of the bezier curve to a monomial polynomial .

		MatrixXd Q, Q_l, Q_u; // Cost Hessian matrix in each block of the objective. No scale, only meta elements .
		MatrixXd MQM; 	      // M' * Q * M in each block of the objective. No scale, only meta elements .
		
		VectorXd C;   // Position coefficients vector, used to record all the pre-compute 'n choose k' combinatorial for the bernstein coefficients .
		VectorXd C_v; // Velocity coefficients vector.
		VectorXd C_a; // Acceleration coefficients vector.
		VectorXd C_j; // Acceleration coefficients vector.

		int poly_num1D = order + 1; 
		M.resize(order + 1, order + 1);
		Q_l = MatrixXd::Zero(order + 1, order + 1);
		Q_u = MatrixXd::Zero(order + 1, order + 1);
		Q   = MatrixXd::Zero(order + 1, order + 1);
		
		C.resize  (order + 1);
		C_v.resize(order    );
		
		if(order > 1)
			C_a.resize(order - 1);

		if(order > 2)
			C_j.resize(order - 2);

		int min_order_l = floor(_min_order);
		int min_order_u = ceil (_min_order);

		if( poly_num1D > min_order_l )
		{
			for( int i = min_order_l; i < poly_num1D; i ++ )
			{
		    	for(int j = min_order_l; j < poly_num1D; j ++ )
		    	{
		      		int coeff = 1.0;
		            int _d = min_order_l - 1;
		                
		            while(_d >= 0)
		            {
		                coeff = coeff * (i - _d) * (j - _d);
		                _d -= 1;
		            }

		            Q_l(i,j) = double(coeff) / double(i + j - 2 * min_order_l + 1);
		    	}
		  	}
		}

		if( poly_num1D > min_order_u )
		{
			for( int i = min_order_u; i < poly_num1D; i ++ )
			{
		    	for(int j = min_order_u; j < poly_num1D; j ++ )
		    	{
		      		int coeff = 1.0;
		            int _d = min_order_u - 1;
		                
		            while(_d >= 0)
		            {
		                coeff = coeff * (i - _d) * (j - _d);
		                _d -= 1;
		            }

		            Q_u(i,j) = double(coeff) / double(i + j - 2 * min_order_u + 1);
		    	}
		  	}
		}

		if(min_order_l == min_order_u)
			Q = Q_u;
		else
			Q = (_min_order - min_order_l) * Q_u + (min_order_u - _min_order) * Q_l;

	 	switch(order)
		{	
			case 0: 
			{
				M << 1;

				break;

			}
			case 1: 
			{
				M << -1,  0,
					 -1,  1;
				break;

			}
			case 2:
			{
				M << -1,  0,  0,
					 -2,  2,  0,
					  1, -2,  1;
				break;

			}
			case 3: 
			{
				M << -1,  0,  0,  0,
					 -3,  3,  0,  0,
					  3, -6,  3,  0,
					 -1,  3, -3,  1;	
				break;

			}
			case 4:
			{
				M <<  1,   0,   0,   0,  0,
					 -4,   4,   0,   0,  0,
					  6, -12,   6,   0,  0,
					 -4,  12, -12,   4,  0,
					  1,  -4,   6,  -4,  1;
				break;
			}
			case 5:
			{
				M << 1,   0,   0,   0,  0,  0,
					-5,   5,   0,   0,  0,  0,
					10, -20,  10,   0,  0,  0,
				   -10,  30, -30,  10,  0,  0,
				     5, -20,  30, -20,  5,  0,
				    -1,   5, -10,  10, -5,  1;
				break;
			}
			case 6:
			{	

				M << 1,   0,   0,   0,   0,  0,  0,
					-6,   6,   0,   0,   0,  0,  0,
					15, -30,  15,   0,   0,  0,  0,
				   -20,  60, -60,  20,   0,  0,  0,
				    15, -60,  90, -60,  15,  0,  0,
				    -6,  30, -60,  60, -30,  6,  0,
				     1,  -6,  15, -20,  15, -6,  1;
				break;
			}
			case 7:
			{
				M << 1,    0,    0,    0,    0,   0,   0,   0,
				    -7,    7,    0,    0,    0,   0,   0,   0,
				    21,   42,   21,    0,    0,   0,   0,   0,
				   -35,  105, -105,   35,    0,   0,   0,   0, 
				    35, -140,  210, -140,   35,   0,   0,   0,
				   -21,  105, -210,  210, -105,  21,   0,   0,
				     7,  -42,  105, -140,  105, -42,   7,   0,
				    -1,    7,  -21,   35,  -35,  21,  -7,   1;
				break;
			}
			case 8:
			{
				M << 1,    0,    0,    0,    0,    0,   0,   0,   0,
				    -8,    8,    0,    0,    0,    0,   0,   0,   0,
				    28,  -56,   28,    0,    0,    0,   0,   0,   0,
				   -56,  168, -168,   56,    0,    0,   0,   0,   0, 
				    70, -280,  420, -280,   70,    0,   0,   0,   0,
				   -56,  280, -560,  560, -280,   56,   0,   0,   0,
				    28, -168,  420, -560,  420, -168,  28,   0,   0,
				    -8,   56, -168,  280, -280,  168, -56,   8,   0,
				     1,   -8,   28,  -56,   70,  -56,  28,  -8,   1;
				break;
			}
			case 9:
			{
				M << 1,    0,     0,     0,     0,    0,    0,     0,     0,    0,
				    -9,    9,     0,     0,     0,    0,    0,     0,     0,    0, 
				    36,  -72,    36,     0,     0,    0,    0,     0,     0,    0, 
				   -84,  252,  -252,    84,     0,    0,    0,     0,     0,    0, 
				   126, -504,   756,  -504,   126,    0,    0,     0,     0,    0,
				  -126,  630, -1260,  1260,  -630,  126,    0,     0,     0,    0,
				    84, -504,  1260, -1680,  1260, -504,   84,     0,     0,    0,
				   -36,  252,  -756,  1260, -1260,  756, -252,    36,     0,    0,
				     9,  -72,   252,  -504,   630, -504,  252,   -72,     9,    0,
				    -1,    9,   -36,    84,  -126,  126,  -84,    36,    -9,    1;
				break;
			}
			case 10:
			{
				M <<  1,     0,     0,     0,      0,     0,    0,     0,     0,    0,   0,
				    -10,    10,     0,     0,      0,     0,    0,     0,     0,    0,   0,
				     45,   -90,    45,     0,      0,     0,    0,     0,     0,    0,   0,
				   -120,   360,  -360,   120,      0,     0,    0,     0,     0,    0,   0,
				    210,  -840,  1260,  -840,    210,     0,    0,     0,     0,    0,   0,
				   -252,  1260, -2520,  2520,  -1260,   252,    0,     0,     0,    0,   0,
				    210, -1260,  3150, -4200,   3150, -1260,  210,     0,     0,    0,   0,
				   -120,  840,  -2520,  4200,  -4200,  2520, -840,   120,     0,    0,   0,
				     45, -360,   1260, -2520,   3150, -2520, 1260,  -360,    45,    0,   0,
				    -10,   90,   -360,   840,  -1260,  1260, -840,   360,   -90,   10,   0,
				      1,  -10,     45,  -120,    210,  -252,  210,  -120,    45,  -10,   1;
				break;
			}
			case 11:
			{
				M <<  1,     0,    0,      0,      0,      0,     0,     0,     0,    0,   0,  0,
				    -11,    11,    0,      0,      0,      0,     0,     0,     0,    0,   0,  0,
				     55,  -110,   55,      0,      0,      0,     0,     0,     0,    0,   0,  0,
				   -165,   495, -495,    165,      0,      0,     0,     0,     0,    0,   0,  0,
				    330, -1320, 1980,  -1320,    330,      0,     0,     0,     0,    0,   0,  0,
				   -462,  2310, -4620,  4620,  -2310,    462,     0,     0,     0,    0,   0,  0,
				    462, -2772,  6930, -9240,   6930,  -2772,   462,     0,     0,    0,   0,  0,
				   -330,  2310, -6930, 11550, -11550,   6930, -2310,   330,     0,    0,   0,  0,
				    165, -1320,  4620, -9240,  11550,  -9240,  4620, -1320,   165,    0,   0,  0,
				    -55,   495, -1980,  4620,  -6930,   6930, -4620,  1980,  -495,   55,   0,  0,
				     11,  -110,   495, -1320,   2310,  -2772,  2310, -1320,   495, -110,  11,  0,
				     -1,    11,   -55,   165,   -330,    462,  -462,   330,  -165,   55, -11,  1;
				break;
			}
			case 12:
			{
				M <<  1,     0,     0,      0,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				    -12,    12,     0,      0,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				     66,  -132,    66,      0,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				   -220,   660,  -660,    220,      0,      0,     0,     0,     0,    0,    0,   0,   0,
				    495, -1980,  2970,  -1980,    495,      0,     0,     0,     0,    0,    0,   0,   0, 
				   -792,  3960, -7920,   7920,  -3960,    792,     0,     0,     0,    0,    0,   0,   0,
				    924, -5544, 13860, -18480,  13860,  -5544,   924,     0,     0,    0,    0,   0,   0,
				   -792,  5544,-16632,  27720, -27720,  16632, -5544,   792,     0,    0,    0,   0,   0,
				    495, -3960, 13860, -27720,  34650, -27720, 13860, -3960,   495,    0,    0,   0,   0,
				   -220,  1980, -7920,  18480, -27720,  27720,-18480,  7920, -1980,  220,    0,   0,   0,
				     66,  -660,  2970,  -7920,  13860, -16632, 13860, -7920,  2970, -660,   66,   0,   0,
				    -12,   132,  -660,   1980,  -3960,   5544, -5544,  3960, -1980,  660, -132,  12,   0,
				      1,   -12,    66,   -220,    495,   -792,   924,  -792,   495, -220,   66, -12,   1;
				break;
			}
		}
		
		MList.push_back(M);
		MQM = M.transpose() * Q * M; // Get the cost block after mapping the coefficients

		MatrixXd F  = CholeskyDecomp(Q);
		MatrixXd FM = F * M;

		MQMList.push_back(MQM);
		FMList.push_back(FM);

		int n = order;
		for(int k = 0; k <= n; k ++ )
		{
			C(k)   = combinatorial(n, k);
			
			if( k <= (n - 1) )
				C_v(k) = combinatorial(n - 1, k);
			if( k <= (n - 2) )
				C_a(k) = combinatorial(n - 2, k);
			if( k <= (n - 3) )
				C_j(k) = combinatorial(n - 3, k);
		}

		CList.push_back(C);
		CvList.push_back(C_v);
		CaList.push_back(C_a);
		CjList.push_back(C_j);

	}
	return ret;
};