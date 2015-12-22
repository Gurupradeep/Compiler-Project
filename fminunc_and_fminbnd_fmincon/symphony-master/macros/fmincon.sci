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



function [xopt,fopt,exitflag,output,lambda,gradient,hessian,zl,zu] = fmincon (varargin)
  // Solves a Constrainted Optimization Problem
  //
  //   Calling Sequence
  //   xopt = fmincon(_f,x0,A,b)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fg)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fg,_fh)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fh)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_cg) 
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fg,_cg)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fg,_fh,_cg)
  //   xopt = fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fh,_cg)
  //   [xopt,fopt] = fmincon(.....)
  //   [xopt,fopt,exitflag]= fmincon(.....)
  //   [xopt,fopt,exitflag,output]= fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient,hessian]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient,hessian,zl]=fmincon(.....)
  //   [xopt,fopt,exitflag,output,lambda,gradient,hessian,zl,zu]=fmincon(.....)
  //
  //   Parameters
  //   _f : a function, representing objective function of the problem 
  //   x0 : a vector of doubles, containing starting of variables of size (1 X n) or (n X 1) where 'n' is the no. of Variables
  //   A : a matrix of doubles, containing coefficients of Linear Inequality Constraints of size (m X n) where 'm' is the no. of Linear Inequality Constraint & 'n' is the 
  //       no. of Variables
  //   b : a vector of doubles, related to 'A' and containing the rhs of the Linear Inequality Constraints of size (m X 1)
  //   Aeq : a matrix of doubles, containing coefficients of Linear Equality Constraints of size (m1 X n) where 'm1' is the no. of Linear Equality Constraint & 'n' is the 
  //         no. of Variables
  //   beq : a vector of doubles, related to 'Aeq' and containing the rhs of the Linear Equality Constraints of size (m1 X 1)
  //   lb : a vector of doubles, containing lower bounds of the variables of size (1 X n) or (n X 1) where 'n' is the no. of Variables
  //   ub : a vector of doubles, containing upperss bounds of the variables of size (1 X n) or (n X 1) where 'n' is the no. of Variables
  //   no_nlic : a scalar of double, related to '_nlc' containing the number of Non-linear In equality constraints in the  Non-linear constraint function('_nlc')
  //   _nlc : a function, representing Non-linear constraint functions(Both Equality and Inequality) of the problem. It is declared in such a way that non-linear 
  //          Inequality constraints are defined first, followed by non-linear Equality constraints
  //		  Note: Constraints should be declared as a vector form (Refer Example Below)  
  //   options: a list, containing option for user to specify -Maximum iteration, Maximum CPU-time, GradObj, HessObj& GradCon.
  //          Syntax for options- options= list("MaxIter", [---], "CpuTime", [---], "GradObj", "ON/OFF", "HessObj", "ON/OFF", "GradCon", "ON/OFF");
  //   		  Default Values for Options==> ("MaxIter", [1000000], "CpuTime", [60], "GradObj", "OFF", "HessObj", "OFF", "GradCon", "OFF");
  //   _fg : a function, representing gradient function of the Objective in Vector Form
  //   _fh : a function, representing hessian function of the Objective in Symmetric Matrix Form
  //   _cg : a function, representing gradient function of the Non-Linear Constraints in Vector Form
  //	     Note: Each element of the constraint's Gradient should be declared seperately as an element of one Vector (Refer Example Below)  
  //   xopt : a vector of doubles, the computed solution of the optimization problem.
  //   fopt : a double, the function value at x
  //   exitflag : Integer identifying the reason the algorithm terminated
  //   output : Structure containing information about the optimization
  //   lambda : a vector of doubles, containing Lagrange multipliers at the optimized point
  //   zl : a vector of doubles, containing lower bound multipliers
  //   zu : a vector of doubles, containing upper bound multipliers
  //   gradient : a vector of doubles, containing Objective's gradient of the optimized point
  //   hessian  : a matrix of doubles, containing Objective's hessian of the optimized point
  //
  //   Description
  //   Search the minimum of a unconstrained optimization problem specified by :
  //   find the minimum of f(x) such that 
  //
  //   <latex>
  //    \begin{eqnarray}
  //    &\mbox{min}_{x}
  //    & f(x) \\
  //    & \text{subject to} & A*x \leq b \\
  //    & & Aeq*x \= beq\\
  //	& & _nlc(x) \leq / = 0\\
  //    & & lb \leq x \leq ub \\
  //    \end{eqnarray}
  //   </latex>
  //
  //   We are calling IPOpt for solving the unconstrained problem, IPOpt is a library written in C++. The code has been written by ​Andreas Wächter and ​Carl Laird.
  //
  //
  // 	Examples
  //	
  //    Example 1:
  //	//Objective function to be minimised
  //	function y=_f(x)
  //		y=x(1)*x(2)+x(2)*x(3);
  //  	endfunction
  //	
  //	//Starting point, linear constraints and variable bounds
  //	x0=[0.1,0.1,0.1];
  //	A=[];
  //	b=[];
  //	Aeq=[];
  //	beq=[];
  //	lb=[]
  //	ub=[]
  //	no_nlic=2;//Number of nonlinear inequality constraints
  //	
  //	//Nonlinear constraints
  //	function [y]=_nlc(x)
  //		y(1)= x(1)^2 - x(2)^2 + x(3)^2 - 2
  //	    	y(2)= x(1)^2 + x(2)^2 + x(3)^2 - 10
  //	
  //	endfunction
  //
  //	//Options
  //	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "ON", "HessObj", "OFF","GradCon", "ON");
  //	
  //	//Gradient of objective function
  //	function y= _funG(x)
  // 		y= [x(2),x(1)+x(3),x(2)];
  //	endfunction
  //	
  //	//Gradient of constraints
  //	function y= _funJ(x)
  //		y(1)= 2*x(1)
  //		y(2)=-2*x(2)
  //		y(3)=2*x(3)
  //		y(4)=2*x(1)
  //		y(5)=2*x(2)
  //		y(6)=2*x(3)
  //		//or, this is another way of writing it : y=[2*x(1),-2*x(2),2*x(3),2*x(1),2*x(2),2*x(3)];	
  //	endfunction
  //	
  //	//Calling fmincon to solve the given problem
  //	[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(_f, x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_funG,_funJ)
  //
  //
  //
  //	Example 2:
  //    //Objective function
  //	function y= _f(x)
  //		y= x(1)^2 + x(2)^2;
  //    endfunction
  //
  //	//Starting point, linear constraints and variable bounds
  //	x0=[1.1,1.1];
  //	A=[];
  //	b=[];
  //	Aeq=[];
  //	beq=[];
  //	lb=[];
  //	ub=[];
  //	no_nlic=1;//Number of nonlinear inequality constraints
  //
  //	//Nonlinear constraints
  //	function y= _nlc(x)
  //		y(1)=x(1)^2 + x(2)^2 -2;//Inequality constraints mentioned before all the equality constraints	     
  //		y(2)= x(1)^2 + x(2)^2 -1;//Equality constraints
  //	endfunction
  //
  //	//Options
  //	options=list("MaxIter", [1000000], "CpuTime", [60], "GradObj", "ON", "HessObj", "OFF", "GradCon", "ON");
  //	//Gradient of objective function
  //	function y=_fg(x)
  //		y=[2*x(1),2*x(2)];
  //	endfunction 
  //
  //	//Gradient of constraints
  //	function [y]=_cg(x)
  //		y=[2*x(1),2*x(2),2*x(1),2*x(2)];	 
  //	endfunction
  //      
  //	//Calling fmincon to solve the given problem
  //	[x,fval,exitflag,output,lambda,grad,hessian]=fmincon(_f,x0,A,b,Aeq,beq,lb,ub,no_nlic,_nlc,options,_fg,_cg)
  //	
  //
  // Authors
  // R.Vidyadhar , Vignesh Kannan
 

	//To check the number of input and output argument
   	[lhs , rhs] = argn();
	
	//To check the number of argument given by user
   	if ( rhs<4 | rhs>14 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,10,11,12,13,14,15"), "fmincon", rhs);
    		error(errmsg)
   	end
    	
	if (rhs==5 | rhs==7 | rhs==9) then
    	errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,10,11,12,13,14,15"), "fmincon", rhs);
    	error(errmsg)
   	end
 
	//Storing the Input Parameters  
   	_f    	 = varargin(1);
   	x0   	 = varargin(2);
   	A    	 = varargin(3);
   	b    	 = varargin(4);
   	Aeq  	 = [];
   	beq  	 = [];
   	lb       = [];
   	ub       = [];
   	no_nlic  =[];
   	_nlc     = [];
   	
   	if (rhs>4) then
   		Aeq  	 = varargin(5);
   		beq  	 = varargin(6);
   	end

   	if (rhs>6) then
   		lb       = varargin(7);
   		ub       = varargin(8);
   	end

   	if (rhs>8) then
   		no_nlic   = varargin(9);
   		_nlc      = varargin(10);
	end
	 
	//To check whether the 1st Input argument (_f) is a function or not
   	if (type(_f) ~= 13 & type(_f) ~= 11) then
   		errmsg = msprintf(gettext("%s: Expected function for Objective (1st Parameter) "), "fmincon");
   		error(errmsg);
   	end
   
	//To check whether the 2nd Input argument (x0) is a Vector/Scalar
   	if (type(x0) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Starting Point (2nd Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check and convert the 2nd Input argument (x0) to row Vector 
   	if((size(x0,1)~=1) & (size(x0,2)~=1)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector or Column Vector for x0 (Starting Point) or Starting Point cannot be Empty"), "fmincon");
   		error(errmsg);
    end

   	if(size(x0,2)==1) then
   		x0=x0';		//Converting x0 to row vector, if it is column vector
   	else 
   	 	x0=x0;		//Retaining the same, if it is already row vector
   	end   	 	
    	s=size(x0);
  	
  	//To check the match between _f (1st Parameter) & x0 (2nd Parameter)
   	if(execstr('init=_f(x0)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "fmincon");
   		error(errmsg);
	end
   	
  	//To check whether the 3rd Input argument (A) is a Matrix/Vector
   	if (type(A) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Matrix/Vector for Constraint Matrix A (3rd parameter)"), "fmincon");
   		error(errmsg);
  	end

	//To check for correct size of A(3rd paramter)
   	if(size(A,2)~=s(2) & size(A,2)~=0) then
   		errmsg = msprintf(gettext("%s: Expected Matrix of size (No of Linear Inequality Constraints X No of Variables) or an Empty Matrix for Linear Inequality Constraint coefficient Matrix A"), "fmincon");
   		error(errmsg);
   	end

   	s1=size(A);
   	
	//To check whether the 4th Input argument (b) is a Vector/Scalar
   	if (type(b) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for b (4th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size of b (4th paramter) and Converting into Column Vector which is required for Ipopt
    if(s1(2)==0) then
    	if(size(b,2)~=0) then
    		errmsg = msprintf(gettext("%s: As Linear Inequality Constraint coefficient Matrix A (3rd parameter) is empty, b (4th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end
	else
   		if((size(b,1)~=1) & (size(b,2)~=1)) then
   			errmsg = msprintf(gettext("%s: Expected Non empty Row/Column Vector for b (4th Parameter) for your Inputs "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)~=s1(1) & size(b,2)==1) then
   			errmsg = msprintf(gettext("%s: Expected Column Vector (No. of Linear inequality Constraints X 1) for b (4th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)==s1(1) & size(b,2)==1) then 
   	 		b=b;
   		elseif(size(b,1)==1 & size(b,2)~=s1(1)) then
   			errmsg = msprintf(gettext("%s: Expected Row Vector (1 X No. of Linear inequality Constraints) for b (4th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)==1 & size(b,2)==s1(1)) then
   			b=b';
   		end 
   	end
  	
  	//To check whether the 5th Input argument (Aeq) is a Matrix/Vector
   	if (type(Aeq) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Matrix/Vector for Equality Constraint Matrix Aeq (5th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size of Aeq (5th paramter)
   	if(size(Aeq,2)~=s(2) & size(Aeq,2)~=0) then
   		errmsg = msprintf(gettext("%s: Expected Matrix of size (No of Linear Equality Constraints X No of Variables) or an Empty Matrix for Linear Equality Constraint coefficient Matrix Aeq"), "fmincon");
   		error(errmsg);
   	end

   	s2=size(Aeq);

	//To check whether the 6th Input argument(beq) is a Vector/Scalar
   	if (type(beq) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for beq (6th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size of beq(6th paramter) and Converting into Column Vector which is required for Ipopt
    if(s2(2)==0) then
    	if(size(beq,2)~=0) then
    		errmsg = msprintf(gettext("%s: As Linear Equality Constraint coefficient Matrix Aeq (5th parameter) is empty, beq (6th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end
   	else
   		if((size(beq,1)~=1) & (size(beq,2)~=1)) then
   			errmsg = msprintf(gettext("%s: Expected Non empty Row/Column Vector for beq (6th Parameter)"), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)~=s2(1) & size(beq,2)==1) then
   			errmsg = msprintf(gettext("%s: Expected Column Vector (No. of Linear Equality Constraints X 1) for beq (6th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)==s2(1) & size(beq,2)==1) then 
   	 		beq=beq;
   		elseif(size(beq,1)==1 & size(beq,2)~=s2(1)) then
   			errmsg = msprintf(gettext("%s: Expected Row Vector (1 X No. of Linear Equality Constraints) for beq (6th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)==1 & size(beq,2)==s2(1)) then
   			beq=beq';
   		end 
   	end
   	
  	
  	//To check whether the 7th Input argument (lb) is a Vector/Scalar
   	if (type(lb) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Lower Bound Vector (7th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check for correct size and data of lb (7th paramter) and Converting it to Column Vector as required by Ipopt
   	if (size(lb,2)==0) then
        	lb = repmat(-%inf,1,s(2));
    end
    
   	if (size(lb,1)~=1) & (size(lb,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound (7th Parameter) should be a vector"), "fmincon");
      error(errmsg); 
    elseif(size(lb,1)~=s(2) & size(lb,2)==1) then
   		errmsg = msprintf(gettext("%s: Expected Column Vector (No. of Variables X 1) for lower bound (7th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(lb,1)==s(2) & size(lb,2)==1) then
   	 	lb=lb;
   	elseif(size(lb,1)==1 & size(lb,2)~=s(2)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector (1 X No. of Variables) for lower bound (7th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(lb,1)==1 & size(lb,2)==s(2)) then
   		lb=lb';
   	end 
   	
   	//To check whether the 8th Input argument (ub) is a Vector/Scalar
   	if (type(ub) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Upper Bound Vector (8th Parameter)"), "fmincon");
   		error(errmsg);
  	end
   	
   	//To check for correct size and data of ub (8th paramter) and Converting it to Column Vector as required by Ipopt
    if (size(ub,2)==0) then
        ub = repmat(%inf,1,s(2));
    end
    
    if (size(ub,1)~=1)& (size(ub,2)~=1) then
      errmsg = msprintf(gettext("%s: Upper Bound (8th Parameter) should be a vector"), "fmincon");
      error(errmsg); 
    elseif(size(ub,1)~=s(2) & size(ub,2)==1) then
   		errmsg = msprintf(gettext("%s: Expected Column Vector (No. of Variables X 1) for upper bound (8th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(ub,1)==s(2) & size(ub,2)==1) then
   	 	ub=ub;
   	elseif(size(ub,1)==1 & size(ub,2)~=s(2)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector (1 X No. of Variables) for beq (8th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(ub,1)==1 & size(ub,2)==s(2)) then
   		ub=ub';
   	end 
    
    //To check the contents of lb & ub (7th & 8th Parameter)
    for i = 1:s(2)
		if (lb(i) == %inf) then
		   	errmsg = msprintf(gettext("%s: Value of Lower Bound can not be infinity"), "fmincon");
    		error(errmsg); 
  		end	

		if (ub(i) == -%inf) then
		   	errmsg = msprintf(gettext("%s: Value of Upper Bound can not be negative infinity"), "fmincon");
    		error(errmsg); 
		end	

		if(ub(i)-lb(i)<=1e-6) then
			errmsg = msprintf(gettext("%s: Difference between Upper Bound and Lower bound should be atleast > 10^6 for variable No.= %d "), "fmincon", i);
    		error(errmsg)
    	end
	end
  	
  	//To check whether the 9th Input argument (no_nlic) is a Scalar
   	if (type(no_nlic) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Scalar for no. of non Linear inequality constraints (9th Parameter)"), "fmincon");
   		error(errmsg);
  	end
  	
  	//To check whether the 10th Input argument (_nlc) is a function or an empty Matrix
   	if (size(no_nlic,1)==1 & size(no_nlic,2)==1) then
   		if (type(_nlc) ~= 13 & type(_nlc) ~= 11) then
   			errmsg = msprintf(gettext("%s: Expected function for non Linear Constraints (10th Parameter) "), "fmincon");
   			error(errmsg);
   		end
   		no_nlc=1;
	elseif (size(no_nlic,2)==0) then	
		no_nlic=0;
		no_nlc=0;
		if (type(_nlc)~=1) then
			errmsg = msprintf(gettext("%s: As no. of Non Linear Inequality Constraints (9th Parameter) is empty, non linear constraint function (10th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end

   		if (size(_nlc,2)~=0) then
   		errmsg = msprintf(gettext("%s: As no. of Non Linear Inequality Constraints (9th Parameter) is empty, non linear constraint function (10th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end
   	else
		errmsg = msprintf(gettext("%s: Expected scalar or Empty matrix for no. of non linear inequality constraints (9th Parameter)"), "fmincon");
    	error(errmsg); 
   	end
	
	//Returns "Invalid Index" Error if size of x0 is not matched with _nlc(10th Parameter)
   	if (no_nlc==1) then
   		if(execstr('init1=_nlc(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Non-Linear Constraint function(9th Parameter) and x0(2nd Parameter) did not match "), "fmincon");
   			error(errmsg);
		end
   		init1=_nlc(x0);
   		no_nlc=size(init1,1);
   		if(no_nlc<no_nlic)
   			errmsg = msprintf(gettext("%s: Error--->Total no. of Non linear Constraints is < than no. of Non linear Inequality Constraints "), "fmincon");
   			error(errmsg);
		end
   	end
   	
   	//To check, Whether Options is been entered by user   
   	if ( rhs<11 ) then
      		param = list();
       else
      		param =varargin(11); //Storing the 3rd Input Parameter in intermediate list named 'param'
    end
   
	//If Options is entered then checking its type for 'list'   
   	if (type(param) ~= 15) then
   		errmsg = msprintf(gettext("%s: Options (11th parameter) should be a list"), "fmincon");
   		error(errmsg);
   	end
   
	//If Options is entered then checking whether even no. of entires are entered   
   	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of Options (list) should be even"), "fmincon");
		error(errmsg);
   	end

	//To set Default Value for Options, If User Doesn't enter Options
   	options = list(..
      		"MaxIter"     , [1000000], ...
      		"CpuTime"   , [60] ...
		);

	//Flags to check whether Gradient is "ON"/"OFF" & Hessian is "ON"/"OFF" 
   	flag1=0;
   	flag2=0;
   	flag3=0;
   	_fg=[];
   	_fh=[];
   	_cg=[];
   	
 
	//To check the User Entry for Options and storing it
   	for i = 1:(size(param))/2
       	select param(2*i-1)
    		case "MaxIter" then
          			options(2*i) = param(2*i);    //Setting the Maximum Iteration as per user entry
       		case "CpuTime" then
          			options(2*i) = param(2*i);    //Setting the Maximum CPU Time as per user entry
        	case "GradObj" then
        				if (param(2*i)=="ON") then
        					//To check whether the user has provided Gradient function if Gradient Option is "ON"
        					if (rhs<12) then      
				     			errmsg = msprintf(gettext("%s: Gradient function of Objective is missing, but GradObj=ON"), "fmincon");
				    			error(errmsg);     			
        					else
        						//This flag1 is activated(ie. =1) if Gradient is supplied
        						flag1=1;
        						pos_fg=12;
        						_fg=varargin(12);
        					end
        				//To check whether Wrong entry(other than ON/OFF) is entered
        				elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        					errmsg = msprintf(gettext("%s: Options for GradObj should be either ON or OFF"), "fmincon");
						error(errmsg);     	
        				end
        	case "HessObj" then
        				if (param(2*i)=="ON") then
        					//To check whether the user has provided Hessian function if Hessian Option is "ON"
							if (flag1==0) then
								if (rhs<12) then    
				     				errmsg = msprintf(gettext("%s: Hessian function of Objective is missing, but HessObj=ON"), "fmincon");
				     				error(errmsg);		
        						else
        							//This flag is activated(ie. =1) if Hessian is supplied
        							flag2=1;
        							pos_fh=12;
        							_fh=varargin(12);
        			    		end         			
        					elseif (flag1==1) then
								if (rhs<13) then    
				     				errmsg = msprintf(gettext("%s: Hessian function of Objective is missing, but HessObj=ON"), "fmincon");
				     				error(errmsg);     			
        						else
        							//This flag is activated(ie. =1) if Hessian is supplied
        							flag2=1;
        							pos_fh=13;
        							_fh=varargin(13);
        			    		end
        			    	end       
        				//To check whether Wrong entry(other than ON/OFF) is entered	            
        				elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        					errmsg = msprintf(gettext("%s: Options for HessObj should be either ON or OFF"), "fmincon");
							error(errmsg);   
        				end	
        	case "GradCon" then
        				if (param(2*i)=="ON") then
        					//To check whether the user has provided Gradient function if Gradient Option is "ON"
        					if (flag1==0 & flag2==0) then
        						if (rhs<12) then      
				     				errmsg = msprintf(gettext("%s: Gradient function of Non-Linear Constraint is missing, but GradCon=ON"), "fmincon");
				    				error(errmsg);     			
        						else
        							pos_cg=12;
        							//This flag is activated(ie. =1) if Hessian is supplied
        							flag3=1;
        							_cg=varargin(12);
        						end
        					elseif((flag1==1 & flag2==0) |(flag1==0 & flag2==1) ) then
        						if (rhs<13) then      
				     				errmsg = msprintf(gettext("%s: Gradient function of Constraints is missing, but GradCon=ON"), "fmincon");
				    				error(errmsg);    			
        						else
        							pos_cg=13;
        							//This flag is activated(ie. =1) if Hessian is supplied
        							flag3=1;
        							_cg=varargin(13);
        						end
        					elseif(flag1==1 & flag2==1) then
        						if (rhs<14) then      
				     				errmsg = msprintf(gettext("%s: Gradient function of Constraints is missing, but GradCon=ON"), "fmincon");
				    				error(errmsg);      			
        						else
        							pos_cg=14;
        							//This flag is activated(ie. =1) if Hessian is supplied
        							flag3=1;
        							_cg=varargin(14);
        						end        					
        					end 	       				      
        			//To check whether Wrong entry(other than ON/OFF) is entered
        				elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        					errmsg = msprintf(gettext("%s: Options for GradCon should be either ON or OFF, but GradCon=ON"), "fmincon");
							error(errmsg);     	
        				end
            else
    	      	errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "fmincon", param(2*i-1));
    	      	error(errmsg)
    		end
   	end
   
   
	//Defining a function to calculate Gradient or Hessian if the respective user entry is OFF 
   	function y=_gradhess(x,t)
		if t==1 then	//To return Gradient
			y=numderivative(_f,x)		
		elseif t==2 then		//To return Hessian]n
			[grad,y]=numderivative(_f,x)
		elseif t==3 then	//To return Gradient
			y=numderivative(_nlc,x)		
		elseif t==4 then		//To return Hessian]n
			[grad,y]=numderivative(_nlc,x)
		end
   	endfunction
 
   //To check the correct no. of inputs given by the user	
   if (flag1==0 & flag2==0 & flag3==0)
   		if(rhs>11) then
        	errmsg = msprintf(gettext("%s: Only 11 Inputs are Needed for this option(GradObj=OFF, HessObj=OFF, GradCon=OFF), but %d were recorded"), "fmincon",rhs);
			error(errmsg); 
		end
   elseif ((flag1==1 & flag2==0 & flag3==0) | (flag1==0 & flag2==1 & flag3==0) | (flag1==0 & flag2==0 & flag3==1)) then
  		if(rhs>12) then
        	errmsg = msprintf(gettext("%s: Only 12 Inputs were needed for this option, but %d were recorded"), "fmincon",rhs);
			error(errmsg);
		end
   elseif ((flag1==1 & flag2==1 & flag3==0) | (flag1==0 & flag2==1 & flag3==1) | (flag1==1 & flag2==0 & flag3==1)) then
   		if(rhs>13) then
        	errmsg = msprintf(gettext("%s: Only 13 Inputs were needed for this option, but %d were recorded"), "fmincon",rhs);
			error(errmsg);
		end
   elseif (flag1==1 & flag2==1 & flag3==1)
   		if(rhs>14) then
        	errmsg = msprintf(gettext("%s: Only 14 Inputs are Needed for this option(GradObj=ON, HessObj=ON, GradCon=ON), but %d were recorded"), "fmincon",rhs);
			error(errmsg); 
		end
	end
	
   //To check the correct input of Gradient and Hessian Functions from Users	     	
   if (flag1==1) then
   		if (type(_fg) ~= 13 & type(_fg) ~= 11) then
  			errmsg = msprintf(gettext("%s: Expected function for Gradient of Objective, since GradObj=ON"), "fmincon");
   			error(errmsg);
   		end
   		if(execstr('sample_fg=_fg(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Objective and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		sample_fg=_fg(x0);
		if (size(sample_fg,1)==s(2) & size(sample_fg,2)==1) then
		elseif (size(sample_fg,1)==1 & size(sample_fg,2)==s(2)) then
		elseif (size(sample_fg,1)~=1 & size(sample_fg,2)~=1) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Gradient function(%dth Parameter)---->Vector function is Expected"), "fmincon",pos_fg);
   			error(errmsg);
   		end
   	end
   	if (flag2==1) then
   		if (type(_fh) ~= 13 & type(_fh) ~= 11) then
  			errmsg = msprintf(gettext("%s: Expected function for Hessian of Objective, since HessObj=ON"), "fmincon");
   			error(errmsg);
   		end
   		if(execstr('sample_fh=_fh(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Hessian function of Objective and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		sample_fh=_fh(x0);
   		if(size(sample_fh,1)~=s(2) | size(sample_fh,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Hessian function(%dth Parameter)---->Symmetric Matrix function is Expected "), "fmincon",pos_fh);
   			error(errmsg);
   		end
   	end
   	if (flag3==1) then
   		if (type(_cg) ~= 13 & type(_cg) ~= 11) then
  			errmsg = msprintf(gettext("%s: Expected function for Gradient of Constraint function,since GradCon=ON"), "fmincon");
   			error(errmsg);
   		end
   		if(execstr('sample_cg=_cg(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Constraint and x0 did not match "), "fmincon", rhs);
   			error(errmsg);
		end
		sample_cg=_cg(x0);
		if (size(sample_cg,1)==no_nlc*s(2) & size(sample_cg,2)==1) then
		elseif (size(sample_cg,1)==1 & size(sample_cg,2)==no_nlc*s(2)) then	
   		elseif(size(sample_cg,1)~=1 & size(sample_cg,2)~=1) then
   			errmsg = msprintf(gettext("%s:  Wrong Input for Constraint Gradient function(%dth Parameter)---->Vector function is Expected (Refer Help)"), "fmincon",pos_cg);
   			error(errmsg);
   		end
   	end
   	
   	//Calling the Ipopt Function for solving the above Problem
    [xopt,fopt,status,iter,cpu,obj_eval,dual,lambda,zl,zu,gradient,hessian1] = solveminconp (_f,_gradhess,A,b,Aeq,beq,lb,ub,no_nlc,no_nlic,_nlc,flag1,_fg,flag2,_fh,flag3,_cg,x0,options)		
   
	//Calculating the values for output   	
   	xopt = xopt';
    	exitflag = status;
    	output = struct("Iterations", [],"Cpu_Time",[],"Objective_Evaluation",[],"Dual_Infeasibility",[]);
   	output.Iterations = iter;
    	output.Cpu_Time = cpu;
    	output.Objective_Evaluation = obj_eval;
    	output.Dual_Infeasibility = dual;

    //Converting hessian of order (1 x (numberOfVariables)^2) received from Ipopt to order (numberOfVariables x numberOfVariables)
    s=size(gradient)
    for i =1:s(2)
    	for j =1:s(2)
			hessian(i,j)= hessian1(j+((i-1)*s(2)))
		end
    end
    
    //In the cases of the problem not being solved return NULL to the output matrices
    if( status~=0 & status~=1 & status~=2 & status~=4 & status~=7 ) then
		xopt=[]
		fopt=[]
		output = struct("Iterations", [],"Cpu_Time",[]);
		output.Iterations = iter;
    		output.Cpu_Time = cpu;
		lambda=[]
		gradient=[]
		hessian=[]
		zl=[]
		zu=[]
    end
		

	//To print Output Message
    select status
    
    	case 0 then
        	printf("\nOptimal Solution Found.\n");
    	case 1 then
        	printf("\nMaximum Number of Iterations Exceeded. Output may not be optimal.\n");
    	case 2 then
       		printf("\nMaximum CPU Time exceeded. Output may not be optimal.\n");
    	case 3 then
        	printf("\nStop at Tiny Step\n");
    	case 4 then
        	printf("\nSolved To Acceptable Level\n");
    	case 5 then
        	printf("\nConverged to a point of local infeasibility.\n");
    	case 6 then
        	printf("\nStopping optimization at current point as requested by user.\n");
    	case 7 then
        	printf("\nFeasible point for square problem found.\n");
    	case 8 then 
        	printf("\nIterates diverging; problem might be unbounded.\n");
    	case 9 then
        	printf("\nRestoration Failed!\n");
    	case 10 then
        	printf("\nError in step computation (regularization becomes too large?)!\n");
    	case 12 then
        	printf("\nProblem has too few degrees of freedom.\n");
    	case 13 then
        	printf("\nInvalid option thrown back by IPOpt\n");
    	case 14 then
        	printf("\nNot enough memory.\n");
    	case 15 then
        	printf("\nINTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
    	else
        	printf("\nInvalid status returned. Notify the Toolbox authors\n");
        	break;
        end
    
    //Remark for the user, If the gradient and hessian is send by the User
    if (no_nlc~=0) then
		disp("||||||Please Make sure you have entered Correct No. of Non-linear Inequality Constraints (9th Parameter) & Non-linear Constraints Functions (10th Parameter) in proper order -->Scilab Will Calculate Based on your input only||||||");	
    end
    
    //Remark for the user, If the gradient and hessian is send by the User
    if (flag1==1 |flag2==1 |flag3==1) then
		disp("||||||Please Make sure you have entered Correct Functions for Gradient or Hessian -->Scilab Will Calculate Based on your input only||||||");
    end
    		
endfunction
