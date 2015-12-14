#include "minuncNLP.hpp"
#include "IpIpoptData.hpp"
#include "sci_iofunc.hpp"

extern "C"
{

#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>
#include <string.h>
#include <assert.h>
#include <iostream>
using namespace std;
//#include <call_scilab.h>

//double x_static,i, *op_obj_x = NULL,*op_obj_value = NULL;

using namespace Ipopt;

minuncNLP::~minuncNLP()
{
	free(finalX_);
	free(finalGradient_);
	free(finalHessian_);
}

//get NLP info such as number of variables,constraints,no.of elements in jacobian and hessian to allocate memory
bool minuncNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	finalGradient_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	finalHessian_ = (double*)malloc(sizeof(double) * numVars_ * numVars_);
	n=numVars_; // Number of variables
	m=numConstr_; // Number of constraints
	nnz_jac_g = 0; // No. of elements in Jacobian of constraints 
	nnz_h_lag = n*(n+1)/2; // No. of elements in lower traingle of Hessian of the Lagrangian.
	index_style=C_STYLE; // Index style of matrices
	return true;
}

//get variable and constraint bound info
bool minuncNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
	unsigned int i;
	for(i=0;i<n;i++)
	{
		x_l[i]=-1.0e19;
		x_u[i]=1.0e19;
	}

        g_l=NULL;
        g_u=NULL;
	return true;
}

// return the value of the constraints: g(x)
bool minuncNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  // return the value of the constraints: g(x)
  g=NULL;
  return true;
}

// return the structure or values of the jacobian
bool minuncNLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
 	if (values == NULL) 
 	{
    		// return the structure of the jacobian of the constraints
    		iRow=NULL; 
    		jCol=NULL;
  	}
  	else 
	{
		values=NULL;
  	}

  	return true;
}

//get value of objective function at vector x
bool minuncNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  	int* funptr=NULL;  
  	if(getFunctionFromScilab(1,&funptr))
  	{
		return 1;
 	}
  	char name[20]="fun_";
  	double obj=0;
  	double *xNew=x;
  	createMatrixOfDouble(pvApiCtx, 3, 1, numVars_, xNew);
  	int positionFirstElementOnStackForScilabFunction = 3;
  	int numberOfRhsOnScilabFunction = 1;
  	int numberOfLhsOnScilabFunction = 1;
  	int pointerOnScilabFunction     = *funptr;
  
  	C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
                               
  	if(getDoubleFromScilab(3,&obj))
  	{
		sciprint("No obj value");
		return 1;
  	}
  	obj_value=obj;  
	
  	return true;
}

//get value of gradient of objective function at vector x.
bool minuncNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  	int* gradhesptr=NULL;
  	if(getFunctionFromScilab(2,&gradhesptr))
  	{
		return 1;
  	}  
  	double *xNew=x;
  	double t=1;
  	createMatrixOfDouble(pvApiCtx, 3, 1, numVars_, xNew);
  	createScalarDouble(pvApiCtx, 4,t);
  	int positionFirstElementOnStackForScilabFunction = 3;
  	int numberOfRhsOnScilabFunction = 2;
  	int numberOfLhsOnScilabFunction = 1;
  	int pointerOnScilabFunction     = *gradhesptr;
	char name[20]="gradhess_";
 
  	C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
 
                               
  	double* resg;  
  	int x0_rows,x0_cols;                           
  	if(getDoubleMatrixFromScilab(3, &x0_rows, &x0_cols, &resg))
  	{
		sciprint("No results");
		return 1;
		
  	}

  	Index i;
  	for(i=0;i<numVars_;i++)
  	{
		grad_f[i]=resg[i];
        	finalGradient_[i]=resg[i];
  	}

  	return true;
}

// This method sets initial values for required vectors . For now we are assuming 0 to all values. 
bool minuncNLP::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
 	assert(init_x == true);
  	assert(init_z == false);
  	assert(init_lambda == false);
	if (init_x == true)
	{ //we need to set initial values for vector x
		for (Index var=0;var<n;var++)
			x[var]=varGuess_[var];//initialize with 0 or we can change.
	}

	return true;
}

/*
 * Return either the sparsity structure of the Hessian of the Lagrangian, 
 * or the values of the Hessian of the Lagrangian  for the given values for
 * x,lambda,obj_factor.
*/
bool minuncNLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
	int* gradhesptr=NULL;
	if(getFunctionFromScilab(2,&gradhesptr))
	{
		return 1;
	}   

	if (values==NULL)
	{
		Index idx=0;
		for (Index row = 0; row < numVars_; row++) 
		{
			for (Index col = 0; col <= row; col++)
			{
				iRow[idx] = row;
				jCol[idx] = col;
				idx++;
		  	}
		}
	}
	else 
	{

		double *xNew=x;
  		double t=2;
	
  		createMatrixOfDouble(pvApiCtx, 3, 1, numVars_, xNew);
  		createScalarDouble(pvApiCtx, 4,t);
  		int positionFirstElementOnStackForScilabFunction = 3;
  		int numberOfRhsOnScilabFunction = 2;
  		int numberOfLhsOnScilabFunction = 1;
  		int pointerOnScilabFunction     = *gradhesptr;
		char name[20]="gradhess_";
  
  		C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
                               
  		double* resh;  
  		int x0_rows,x0_cols;                           
  		if(getDoubleMatrixFromScilab(3, &x0_rows, &x0_cols, &resh))
		{
			sciprint("No results");
			return 1;
		}

		Index index=0;
		for (Index row=0;row < numVars_ ;++row)
		{
			for (Index col=0; col <= row; ++col)
			{
				values[index++]=obj_factor*(resh[numVars_*row+col]);
			}
		}

		Index i;
  		for(i=0;i<numVars_*numVars_;i++)
  		{
        		finalHessian_[i]=resh[i];
 		}
	}
	
       	return true;
}


void minuncNLP::finalize_solution(SolverReturn status,Index n, const Number* x, const Number* z_L, const Number* z_U,Index m, const Number* g, const Number* lambda, Number obj_value,const IpoptData* ip_data,IpoptCalculatedQuantities* ip_cq)
{
	finalX_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	for (Index i=0; i<numVars_; i++) 
	{
    		 finalX_[i] = x[i];
	}

	finalObjVal_ = obj_value;
	status_ = status;
	if (status_ == 0 | status_ == 1 | status_ == 2)
	{
		iter_ = ip_data->iter_count();
	}
}


const double * minuncNLP::getX()
{	
	return finalX_;
}

const double * minuncNLP::getGrad()
{	
	return finalGradient_;
}

const double * minuncNLP::getHess()
{	
	return finalHessian_;
}

double minuncNLP::getObjVal()
{	
	return finalObjVal_;
}

double minuncNLP::iterCount()
{	
	return (double)iter_;
}

int minuncNLP::returnStatus()
{	
	return status_;
}

}

