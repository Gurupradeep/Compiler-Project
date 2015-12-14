//NAME    : VIGNESH K
//DATE    : 18/11/2015
//COLLEGE : NATIONAL INSTITUTE OF TECHNOLOGY KARNATAKA
//BRANCH  : COMPUTER SCIENCE AND ENGINEERING, B.Tech 2nd YEAR 

#ifndef __minuncNLP_HPP__
#define __minuncNLP_HPP__
#include "IpTNLP.hpp"

using namespace Ipopt;

class minuncNLP : public TNLP
{
	private:

  	Index numVars_;	                 //Number of input variables

  	Index numConstr_;                //Number of constraints 

  	const Number *varGuess_= NULL;	 //varGuess_ is a pointer to a matrix of size of 1*numVars_
				         //with initial guess of all variables.

  	Number *finalX_= NULL;           //finalX_ is a pointer to a matrix of size of 1*numVars_
				         //with final value for the primal variables.

  	Number *finalGradient_=NULL;     //finalGradient_ is a pointer to a matrix of size of numVars_*numVars_
				         //with final value of gradient for the primal variables.


  	Number *finalHessian_=NULL;      //finalHessian_ is a pointer to a matrix of size of 1*numVar_
				         //with final value of hessian for the primal variables.


  	Number finalObjVal_;          	 //finalObjVal_ is a scalar with the final value of the objective.

  	int iter_;			 //Number of iteration.

  	int status_;			 //Solver return status


  	minuncNLP(const minuncNLP&);
  	minuncNLP& operator=(const minuncNLP&);

	public:

  	/** user defined constructor */
  	minuncNLP(Index nV, Index nC,Number *x0):numVars_(nV),numConstr_(nC),varGuess_(x0),finalX_(0),finalGradient_(0),finalHessian_(0),finalObjVal_(1e20){	}

  	/** default destructor */
  	virtual ~minuncNLP();

  	/** Method to return some info about the nlp */
  	virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  	/** Method to return the bounds for my problem */
  	virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  	/** Method to return the starting point for the algorithm */
  	virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  	/** Method to return the objective value */
  	virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  	/** Method to return the gradient of the objective */
  	virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  	/** Method to return the constraint residuals */
  	virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  	/** Method to return:
  	*   1) The structure of the jacobian (if "values" is NULL)
   	*   2) The values of the jacobian (if "values" is not NULL)
   	*/
  	virtual bool eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values);

  	/** Method to return:
   	*   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   	*   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   	*/
  	virtual bool eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values);

  	/** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  	virtual void finalize_solution(SolverReturn status,Index n, const Number* x, const Number* z_L, const Number* z_U,Index m, const Number* g, const Number* lambda,Number obj_value,const IpoptData* ip_data,IpoptCalculatedQuantities* ip_cq);
  
  	const double * getX();		//Returns a pointer to a matrix of size of 1*numVars_ 
					//with final value for the primal variables.
  
  	const double * getGrad();       //Returns a pointer to a matrix of size of 1*numVars_ 
					//with final value of gradient for the primal variables.

  	const double * getHess();       //Returns a pointer to a matrix of size of numVars_*numVars_ 
					//with final value of hessian for the primal variables.

  	double getObjVal();		//Returns the output of the final value of the objective.

  	double iterCount();		//Returns the iteration count

  	int returnStatus();		//Returns the status count

};


#endif
