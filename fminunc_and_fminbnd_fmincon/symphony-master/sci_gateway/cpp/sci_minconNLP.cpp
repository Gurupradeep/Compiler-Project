// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: R.Vidyadhar & Vignesh Kannan
// Organization: FOSSEE, IIT Bombay
// Email: rvidhyadar@gmail.com & vignesh2496@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


#include "minconNLP.hpp"
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
using namespace Ipopt;

minconNLP::~minconNLP()
{
	free(finalX_);
	free(finalGradient_);
	free(finalHessian_);
	free(finalZu_);
	free(finalZl_);
	free(finalLambda_);
}

//get NLP info such as number of variables,constraints,no.of elements in jacobian and hessian to allocate memory
bool minconNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	finalGradient_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	finalHessian_ = (double*)malloc(sizeof(double) * numVars_ * numVars_);
	
	n=numVars_; // Number of variables
	m=numConstr_; // Number of constraints

	nnz_jac_g = n*m; // No. of elements in Jacobian of constraints 
	nnz_h_lag = n*(n+1)/2; // No. of elements in lower traingle of Hessian of the Lagrangian.

	index_style=C_STYLE; // Index style of matrices
	return true;
}

//get variable and constraint bound info
bool minconNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
	unsigned int i;
	
	//assigning bounds for the variables
	for(i=0;i<n;i++)
	{
		x_l[i]=varLB_[i];
		x_u[i]=varUB_[i];
	}
	
	if(m==0)
        {
		g_l=NULL;
	        g_u=NULL;
	}

	else
	{
		unsigned int c=0;

		//bounds of linear equality constraints
		for(i=0;i<Aeqrows_;i++)
		{
			g_l[c]=g_u[c]=beq_[i];
			c++;
		}

		//bounds of linear inequality constraints
		for(i=0;i<Arows_;i++)
		{
			g_l[c]=-1.0e19;
			g_u[c]=b_[i];
			
			c++;
		}

		//bounds of non-linear inequality constraints
		for(i=0;i<nonlinIneqCon_;i++)
		{
			g_l[c]=-1.0e19;
			g_u[c]=0;
			c++;
		}
	
		//bounds of non-linear equality constraints
		for(i=0;i<nonlinCon_-nonlinIneqCon_;i++)
		{
			g_l[c]=g_u[c]=0;
			c++;
		}
	}

	return true;
}

// return the value of the constraints: g(x)
bool minconNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  	// return the value of the constraints: g(x)

	unsigned int i;
	unsigned int j;

	if(m==0)
  		g=NULL;

	else
	{
		unsigned int c=0;

		//value of linear equality constraints
		for(i=0;i<Aeqrows_;i++)
		{
			g[c]=0;
			for(j=0;j<Aeqcols_;j++)
				g[c] += Aeq_[j*Aeqrows_+i]*x[j];
			c++;			
		}

		//value of linear inequality constraints
		for(i=0;i<Arows_;i++)
		{
			g[c]=0;
			for(j=0;j<Acols_;j++)
				g[c] += A_[j*Arows_+i]*x[j];
			c++;
		}

		//value of non-linear constraints
		if(nonlinCon_!=0)
		{
			int* constr=NULL;  
	  		if(getFunctionFromScilab(11,&constr))
	  		{
				return 1;
	 		}
	  		char name[20]="_nlc";
		  	double *xNew=x;
  			createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  			int positionFirstElementOnStackForScilabFunction = 18;
  			int numberOfRhsOnScilabFunction = 1;
  			int numberOfLhsOnScilabFunction = 1;
  			int pointerOnScilabFunction     = *constr;
  
 		 	C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));

			double* resc;  
		  	int xC_rows,xC_cols;                           
		  	if(getDoubleMatrixFromScilab(18, &xC_rows, &xC_cols, &resc))
		  	{
				sciprint("No results");
				return 1;
		
		  	}

			for(i=0;i<nonlinCon_;i++)
			{
				g[c]=resc[i];
				c++;
			}
		}
	}

  	return true;
}

// return the structure or values of the jacobian
bool minconNLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
 	if (values == NULL) 
 	{
    		if(m==0)// return the structure of the jacobian of the constraints
    		{
			iRow=NULL; 
    			jCol=NULL;
  		}

		else
		{
			unsigned int i,j,idx=0;
			for(int i=0;i<m;i++)
				for(j=0;j<n;j++)
				{
					iRow[idx]=i;
					jCol[idx]=j;
					idx++;
				}
		}
	}
					
  	else 
	{
		if(m==0)
			values=NULL;

		else
		{
			unsigned int i,j,c=0;

			//jacobian of linear equality constraints
			for(i=0;i<Aeqrows_;i++)
			{
				for(j=0;j<Aeqcols_;j++)
				{
					values[c] = Aeq_[j*Aeqrows_+i];
					c++;
				}
			}

			//jacobian of linear inequality constraints
			for(i=0;i<Arows_;i++)
			{
				for(j=0;j<Acols_;j++)
				{
					values[c] = A_[j*Arows_+i];
					c++;
				}
			}

			//jacobian of non-linear constraints
			if(nonlinCon_!=0)
			{
				if(flag3_==0)
				{
					int* gradhessptr=NULL;
  					if(getFunctionFromScilab(2,&gradhessptr))
  					{
						return 1;
  					}  
  					double *xNew=x;
  					double t=3;
  					createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  					createScalarDouble(pvApiCtx, 19,t);
  					int positionFirstElementOnStackForScilabFunction = 18;
  					int numberOfRhsOnScilabFunction = 2;
  					int numberOfLhsOnScilabFunction = 1;
  					int pointerOnScilabFunction     = *gradhessptr;
					char name[20]="_gradhess";
 	
  					C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                        	                                       &numberOfLhsOnScilabFunction,
                        	                                       &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
	
					double* resj;  
  					int xJ_rows,xJ_cols;                           
  					if(getDoubleMatrixFromScilab(18, &xJ_rows, &xJ_cols, &resj))
  					{
						sciprint("No results");
						return 1;
					}
	
				        for(i=0;i<nonlinCon_;i++)
					{
						for(j=0;j<n;j++)
						{
							values[c] = resj[j*(int)nonlinCon_+i];
							c++;
						}
					}			
				}
			
				else
				{
					int* jacptr=NULL;
  					if(getFunctionFromScilab(17,&jacptr))
  					{
						return 1;
  					}  
	
  					double *xNew=x;
  					createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  					int positionFirstElementOnStackForScilabFunction = 18;
  					int numberOfRhsOnScilabFunction = 1;
  					int numberOfLhsOnScilabFunction = 1;
  					int pointerOnScilabFunction     = *jacptr;
					char name[20]="_cg";
 	
  					C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                        	                                       &numberOfLhsOnScilabFunction,
                        	                                       &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
	
					double* resj;  
  					int xJ_rows,xJ_cols;                           
  					if(getDoubleMatrixFromScilab(18, &xJ_rows, &xJ_cols, &resj))
  					{
						sciprint("No results");
						return 1;
					}
	
				        for(i=0;i<nonlinCon_*n;i++)
					{
						values[c] = resj[i];
						c++;
					}			
				}
			}
		}	
	}	

  	return true;
}

//get value of objective function at vector x
bool minconNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  	int* funptr=NULL;  
  	if(getFunctionFromScilab(1,&funptr))
  	{
		return 1;
 	}
  	char name[20]="_f";
  	double obj=0;
  	double *xNew=x;
  	createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  	int positionFirstElementOnStackForScilabFunction = 18;
  	int numberOfRhsOnScilabFunction = 1;
  	int numberOfLhsOnScilabFunction = 1;
  	int pointerOnScilabFunction     = *funptr;
  
  	C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
                               
  	if(getDoubleFromScilab(18,&obj))
  	{
		sciprint("No obj value");
		return 1;
  	}
  	obj_value=obj;  
	
  	return true;
}

//get value of gradient of objective function at vector x.
bool minconNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  	if (flag1_==0)
  	{	
  		int* gradhessptr=NULL;
  		if(getFunctionFromScilab(2,&gradhessptr))
  		{
			return 1;
  		}  
  		double *xNew=x;
  		double t=1;
  		createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  		createScalarDouble(pvApiCtx, 19,t);
  		int positionFirstElementOnStackForScilabFunction = 18;
  		int numberOfRhsOnScilabFunction = 2;
  		int numberOfLhsOnScilabFunction = 1;
  		int pointerOnScilabFunction     = *gradhessptr;
		char name[20]="_gradhess";
 
  		C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
  	}

  	else 
  	{
  		int* gradptr=NULL;
  		if(getFunctionFromScilab(13,&gradptr))
  		{
			return 1;
  		}  
  		double *xNew=x;
	  	createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  		int positionFirstElementOnStackForScilabFunction = 18;
  		int numberOfRhsOnScilabFunction = 1;
  		int numberOfLhsOnScilabFunction = 1;
  		int pointerOnScilabFunction     = *gradptr;
		char name[20]="_fg";
 	
  		C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
        	                                                       &numberOfLhsOnScilabFunction,
        	                                                       &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
   	}

	double* resg;  
  	int x0_rows,x0_cols;                           
  	if(getDoubleMatrixFromScilab(18, &x0_rows, &x0_cols, &resg))
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
bool minconNLP::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
 	assert(init_x == true);
  	assert(init_z == false);
  	assert(init_lambda == false);
	if (init_x == true)
	{ //we need to set initial values for vector x
		for (Index var=0;var<n;var++)
			x[var]=varGuess_[var];
	}

	return true;
}

/*
 * Return either the sparsity structure of the Hessian of the Lagrangian, 
 * or the values of the Hessian of the Lagrangian  for the given values for
 * x,lambda,obj_factor.
*/

bool minconNLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
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
		//hessian of the objective function		
		if(flag2_==0)
	  	{
			int* gradhessptr=NULL;
			if(getFunctionFromScilab(2,&gradhessptr))
			{
				return 1;
			}          	
			double *xNew=x;
  			double t=2;
			createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  			createScalarDouble(pvApiCtx, 19,t);
  			int positionFirstElementOnStackForScilabFunction = 18;
  			int numberOfRhsOnScilabFunction = 2;
  			int numberOfLhsOnScilabFunction = 1;
  			int pointerOnScilabFunction     = *gradhessptr;
			char name[20]="_gradhess";
  	
  			C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                	                                               &numberOfLhsOnScilabFunction,
                	                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
 	    	}	

 	    	else
 	    	{		
			int* hessptr=NULL;
			if(getFunctionFromScilab(15,&hessptr))
			{
				return 1;
			}          	
			double *xNew=x;	
  			createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
  			int positionFirstElementOnStackForScilabFunction = 18;
  			int numberOfRhsOnScilabFunction = 1;
  			int numberOfLhsOnScilabFunction = 1;
  			int pointerOnScilabFunction     = *hessptr;
			char name[20]="_fh";
  	
  			C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                	                                               &numberOfLhsOnScilabFunction,
                	                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));	
 	    	}	

	        double* resTemph;  
  		int x0_rows,x0_cols;                           
  		if(getDoubleMatrixFromScilab(18, &x0_rows, &x0_cols, &resTemph))
		{
			sciprint("No results");
			return 1;
		}

		double* resh=(double*)malloc(sizeof(double)*n*n);
		Index i;
		for(i=0;i<numVars_*numVars_;i++)
		{
       			resh[i]=resTemph[i];
		}

		//sum of hessians of constraints each multiplied by its own lambda factor
		double* sum=(double*)malloc(sizeof(double)*n*n);
		if(nonlinCon_!=0)
		{	
					
			int* gradhessptr=NULL;
			if(getFunctionFromScilab(2,&gradhessptr))
			{
				return 1;
			}          	
			
			double *xNew=x;
	  		double t=4;
			createMatrixOfDouble(pvApiCtx, 18, 1, numVars_, xNew);
	  		createScalarDouble(pvApiCtx, 19,t);
	  		int positionFirstElementOnStackForScilabFunction = 18;
	  		int numberOfRhsOnScilabFunction = 2;
	  		int numberOfLhsOnScilabFunction = 1;
	  		int pointerOnScilabFunction     = *gradhessptr;
			char name[20]="_gradhess";
	  	
	  		C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
	                	                                               &numberOfLhsOnScilabFunction,
	                	                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
	
	  		double* resCh;
	  		int xCh_rows,xCh_cols;                           
	  		if(getDoubleMatrixFromScilab(18, &xCh_rows, &xCh_cols, &resCh))
			{
				sciprint("No results");
				return 1;
			}
		
			Index j;
			
			for(i=0;i<numVars_*numVars_;i++)
			{
				sum[i]=0;
				for(j=0;j<nonlinCon_;j++)
					sum[i]+=lambda[j+Aeqrows_+Arows_]*resCh[i*(int)nonlinCon_+j];
			}
		}
		
		else
		{
			for(i=0;i<numVars_*numVars_;i++)	
				sum[i]=0;
		}
		
		//computing the lagrangian
		Index index=0;
		for (Index row=0;row < numVars_ ;++row)
		{
			for (Index col=0; col <= row; ++col)
			{
				values[index++]=obj_factor*(resh[numVars_*row+col])+sum[numVars_*row+col];
			}
		}

		for(i=0;i<numVars_*numVars_;i++)
		{
       			finalHessian_[i]=resh[i];
		}	

		free(resh);
		free(sum);
	}	

       	return true;
}

//returning the results
void minconNLP::finalize_solution(SolverReturn status,Index n, const Number* x, const Number* z_L, const Number* z_U,Index m, const Number* g, const Number* lambda, Number obj_value,const IpoptData* ip_data,IpoptCalculatedQuantities* ip_cq)
{
	finalX_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	for (Index i=0; i<numVars_; i++) 
	{
    		 finalX_[i] = x[i];
	}

	finalZl_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	for (Index i=0; i<n; i++) 
	{
    		 finalZl_[i] = z_L[i];
	}

	finalZu_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	for (Index i=0; i<n; i++) 
	{
    		 finalZu_[i] = z_U[i];
	}

	finalLambda_ = (double*)malloc(sizeof(double) * numConstr_ * 1);
	for (Index i=0; i<m; i++) 
	{
    		 finalLambda_[i] = lambda[i];
	}

	finalObjVal_ = obj_value;
	status_ = status;
	iter_ = ip_data->iter_count();
}


const double * minconNLP::getX()
{	
	return finalX_;
}

const double * minconNLP::getGrad()
{	
	return finalGradient_;
}

const double * minconNLP::getHess()
{	
	return finalHessian_;
}

const double * minconNLP::getZl()
{	
	return finalZl_;
}

const double * minconNLP::getZu()
{	
	return finalZu_;
}

const double * minconNLP::getLambda()
{	
	return finalLambda_;
}

double minconNLP::getObjVal()
{	
	return finalObjVal_;
}

double minconNLP::iterCount()
{	
	return (double)iter_;
}

int minconNLP::returnStatus()
{	
	return status_;
}

}



