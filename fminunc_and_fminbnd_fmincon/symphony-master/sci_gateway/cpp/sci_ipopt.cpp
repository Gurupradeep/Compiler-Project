/*
 * Quadratic Programming Toolbox for Scilab using IPOPT library
 * Authors :
	Sai Kiran
	Keyur Joshi
	Iswarya
	Harpreet Singh
 */

#include "sci_iofunc.hpp"
#include "IpIpoptApplication.hpp"
#include "QuadNLP.hpp"

extern "C"{
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>

int sci_solveqp(char *fname)
{
	using namespace Ipopt;

	CheckInputArgument(pvApiCtx, 11, 11); // We need total 11 input arguments.
	CheckOutputArgument(pvApiCtx, 7, 7);
	
	// Error management variable
	SciErr sciErr;

	// Input arguments
	double *QItems=NULL,*PItems=NULL,*ConItems=NULL,*conUB=NULL,*conLB=NULL;
	double *cpu_time=NULL,*max_iter=NULL,*varUB=NULL,*varLB=NULL,*init_guess=NULL;
	static unsigned int nVars = 0,nCons = 0;
	unsigned int temp1 = 0,temp2 = 0, iret = 0;

	// Output arguments
	double *fX = NULL, ObjVal=0,iteration=0, *Zl=NULL, *Zu=NULL, *Lambda=NULL;
	int rstatus = 0;

	////////// Manage the input argument //////////
	
	//Number of Variables
	if(getIntFromScilab(1,&nVars))
	{
		return 1;
	}

	//Number of Constraints
	if (getIntFromScilab(2,&nCons))
	{
		return 1;
	}

	//Q matrix from scilab
	temp1 = nVars;
	temp2 = nCons;
	if (getFixedSizeDoubleMatrixFromScilab(3,temp1,temp1,&QItems))
	{
		return 1;
	}
	
	//P matrix from scilab
	temp1 = 1;
	temp2 = nVars; 
	if (getFixedSizeDoubleMatrixFromScilab(4,temp1,temp2,&PItems))
	{
		return 1;
	}

	if (nCons!=0)
	{
		//conMatrix matrix from scilab
		temp1 = nCons;
		temp2 = nVars;

		if (getFixedSizeDoubleMatrixFromScilab(5,temp1,temp2,&ConItems))
		{
			return 1;
		}

		//conLB matrix from scilab
		temp1 = 1;
		temp2 = nCons;
		if (getFixedSizeDoubleMatrixFromScilab(6,temp1,temp2,&conLB))
		{
			return 1;
		}

		//conUB matrix from scilab
		if (getFixedSizeDoubleMatrixFromScilab(7,temp1,temp2,&conUB))
		{
			return 1;
		}
	}

	//varLB matrix from scilab
	temp1 = 1;
	temp2 = nVars;
	if (getFixedSizeDoubleMatrixFromScilab(8,temp1,temp2,&varLB))
	{
		return 1;
	}


	//varUB matrix from scilab
	if (getFixedSizeDoubleMatrixFromScilab(9,temp1,temp2,&varUB))
	{
		return 1;
	}

	//Initial Value of variables from scilab
	if (getFixedSizeDoubleMatrixFromScilab(10,temp1,temp2,&init_guess))
	{
		return 1;
	}
	
	//Getting the parameters

	temp1 = 1;
	temp2 = 1;

	//Getting maximum iteration
	if (getFixedSizeDoubleMatrixInList(11,2,temp1,temp2,&max_iter))
	{
		return 1;
	}

	//Getting Cpu Time
	if (getFixedSizeDoubleMatrixInList(11,4,temp1,temp2,&cpu_time))
	{
		return 1;
	}

	// Starting Ipopt

	SmartPtr<QuadNLP> Prob = 
	new QuadNLP(nVars,nCons,QItems,PItems,ConItems,conUB,conLB,varUB,varLB,init_guess);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->RethrowNonIpoptException(true);

	////////// Managing the parameters //////////

	app->Options()->SetNumericValue("tol", 1e-7);
	app->Options()->SetIntegerValue("max_iter", (int)*max_iter);
	app->Options()->SetNumericValue("max_cpu_time", *cpu_time);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	// Indicates whether all equality constraints are linear 
	app->Options()->SetStringValue("jac_c_constant", "yes");
	// Indicates whether all inequality constraints are linear 
	app->Options()->SetStringValue("jac_d_constant", "yes");	
	// Indicates whether the problem is a quadratic problem 
	app->Options()->SetStringValue("hessian_constant", "yes");

	///////// Initialize the IpoptApplication and process the options /////////
	ApplicationReturnStatus status;
 	status = app->Initialize();
	if (status != Solve_Succeeded) {
	  	sciprint("\n*** Error during initialization!\n");
   	 return (int) status;
 	 }
	 // Ask Ipopt to solve the problem
	
	 status = app->OptimizeTNLP(Prob);

	rstatus = Prob->returnStatus();

	////////// Manage the output argument //////////

	if (rstatus == 0 | rstatus == 1 | rstatus == 2){
		fX = Prob->getX();
		ObjVal = Prob->getObjVal();
		iteration = Prob->iterCount();

		if (returnDoubleMatrixToScilab(1, 1, nVars, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(4, 1, 1, &iteration))
		{
			return 1;
		}
	}

	else
	{
		if (returnDoubleMatrixToScilab(1, 0, 0, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(4, 1, 1, &iteration))
		{
			return 1;
		}
	}


	if(rstatus == 0){
	
		Zl = Prob->getZl();
		Zu = Prob->getZu();
		Lambda = Prob->getLambda();
		if (returnDoubleMatrixToScilab(5, 1, nVars, Zl))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(6, 1, nVars, Zu))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(7, 1, 1, Lambda))
		{
			return 1;
		}
	}

	else{
		if (returnDoubleMatrixToScilab(5, 0, 0, Zl))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(6, 0, 0, Zu))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(7, 0, 0, Lambda))
		{
			return 1;
		}
	}
	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically
	// be deleted.

	return 0;
	}
}
