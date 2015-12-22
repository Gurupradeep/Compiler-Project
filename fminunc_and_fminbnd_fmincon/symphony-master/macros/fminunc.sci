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



function [xopt,fopt,exitflag,output,gradient,hessian] = fminunc (varargin)
  // Solves a Unconstrainted Optimization Problem
  //
  //   Calling Sequence
  //   xopt = fminunc(_f,x0)
  //   xopt = fminunc(_f,x0,options)
  //   xopt = fminunc(_f,x0,options,_g)
  //   xopt = fminunc(_f,x0,options,_h)
  //   xopt = fminunc(_f,x0,options,_g,_h)
  //   [xopt,fopt] = fminunc(.....)
  //   [xopt,fopt,exitflag]= fminunc(.....)
  //   [xopt,fopt,exitflag,output]= fminunc(.....)
  //   [xopt,fopt,exitflag,output,gradient]=fminunc(.....)
  //   [xopt,fopt,exitflag,output,gradient,hessian]=fminunc(.....)
  //
  //   Parameters
  //   _f : a function, representing objective function of the problem 
  //   x0 : a vector of doubles, containing starting of variables.
  //   options: a list, containing option for user to specify -Maximum iteration, Maximum CPU-time, Gradient&  Hessian
  //            Syntax for options- options= list("MaxIter", [---], "CpuTime", [---], "Gradient", "ON/OFF", "Hessian", "ON/OFF");
  //   		    Default Values for Options==> ("MaxIter", [1000000], "CpuTime", [60], "Gradient", "OFF", "Hessian", "OFF");
  //   _g : a function, representing gradient function of the problem in Vector Form 
  //   _h : a function, representing hessian function of the problem in Symmetric Matrix form
  //   xopt : a vector of doubles, the computed solution of the optimization problem.
  //   fopt : a double, the function value at x.
  //   exitflag : Integer identifying the reason the algorithm terminated.
  //   output : Structure containing information about the optimization.
  //   gradient : a vector of doubles, containing the gradient of the optimized point.
  //   hessian  : a matrix of doubles, containing the hessian of the optimized point.
  //
  //   Description
  //   Search the minimum of a unconstrained optimization problem specified by :
  //   find the minimum of f(x) such that 
  //
  //   <latex>
  //    \begin{eqnarray}
  //    &\mbox{min}_{x}
  //    & f(x)\\
  //    \end{eqnarray}
  //   </latex>
  //
  //   We are calling IPOpt for solving the unconstrained problem, IPOpt is a library written in C++. The code has been written by ​Andreas Wächter and ​Carl Laird.
  //
  //    Examples
  //      //Find x in R^2 such that it minimizes rosenbrock function 
  //      //f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2
  //
  //      function y= _f(x)
  //   	     y= 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
  //      endfunction
  //      function y= _g(x)
  //   	     y= [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1)), 200*(x(2)-x(1)^2)]; //Row Vector is expected for gradient function
  //     endfunction
  //     function y= _h(x)
  //   	     y= [1200*x(1)^2, -400*x(1);-400*x(1), 200 ]; //symmentric Matrix is expected for hessian function
  //     endfunction
  //     x0=[2,7];
  //     options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", "ON", "Hessian", "ON");
  //     [xopt,fopt,exitflag,output,gradient,hessian]=fminunc(_f,x0,options,_g,_h)
  //
  //
  // Examples
  //      //Find x in R^2 such that the below function is minimum
  //      //f = x(1)^2 + x(2)^2
  //
  //      function y= _f(x)
  //   	     y= x(1)^2 + x(2)^2;
  //      endfunction
  //      x0=[2,1];
  //      [xopt,fopt]=fminunc(_f,x0)
  //
  // Authors
  // R.Vidyadhar , Vignesh Kannan
    

	//To check the number of input and output argument
   	[lhs , rhs] = argn();
	
	//To check the number of argument given by user
   	if ( rhs<2 | rhs>5 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2 or 5"), "fminunc", rhs);
    		error(errmsg)
   	end
 
	//Storing the 1st and 2nd Input Parameters  
   	_f = varargin(1);
   	x0 = varargin(2);
      
	//To check whether the 1st Input argument(_f) is function or not
   	if (type(_f) ~= 13 & type(_f) ~= 11) then
   		errmsg = msprintf(gettext("%s: Expected function for Objective "), "fminunc");
   		error(errmsg);
   	end
   
	//To check whether the 2nd Input argument(x0) is Vector/Scalar
   	if (type(x0) ~= 1) then
   		errmsg = msprintf(gettext("%s: Expected Vector/Scalar for Starting Point"), "fminunc");
   		error(errmsg);
  	end
   
	//To check and convert the 2nd Input argument(x0) to row Vector 
   	if((size(x0,1)~=1) & (size(x0,2)~=1)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector or Column Vector for x0 (Initial Value) "), "fminunc", rhs);
   		error(errmsg);
   	else
   		if(size(x0,2)==1) then
   			x0=x0';		//Converting x0 to row vector, if it is column vector
   		else 
   	 		x0=x0;		//Retaining the same, if it is already row vector
   		end   	 	
        	s=size(x0);	
   	end
   

  	//To check the match between _f (1st Parameter) & x0 (2nd Parameter)
   	if(execstr('init=_f(x0)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "fminunc");
   		error(errmsg);
	end
   
	//To check, Whether Options is been entered by user   
   	if ( rhs<3  ) then
      		param = list();
       else
      		param =varargin(3); //Storing the 3rd Input Parameter in intermediate list named 'param'
    
   	end
   
	//If Options is entered then checking its type for 'list'   
   	if (type(param) ~= 15) then
   		errmsg = msprintf(gettext("%s: 3rd Input parameter should be a list (ie. Options) "), "fminunc");
   		error(errmsg);
   	end
   
	//If Options is entered then checking whether even no. of entires are entered   
   	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "fminunc");
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
   	_g=[]
   	_h=[]
 
	//To check the User Entry for Options and storing it
   	for i = 1:(size(param))/2
       	select param(2*i-1)
    		case "MaxIter" then
          			options(2*i) = param(2*i);    //Setting the Maximum Iteration as per user entry
       		case "CpuTime" then
          			options(2*i) = param(2*i);    //Setting the Maximum CPU Time as per user entry
        	case "Gradient" then
        			if (param(2*i)=="ON") then
        				//To check whether the user has provided Gradient function if Gradient Option is "ON"
        				if (rhs<4) then      
				     		errmsg = msprintf(gettext("%s: Gradient function is missing"), "fminunc");
				    		error(errmsg);     			
        				end
        				//This flag is activated(ie. =1) if Gradient is supplied
        				flag1=1;
        				pos_g=4;
        				_g=varargin(4);        				      
        			//To check whether Wrong entry(other than ON/OFF) is entered
        			elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        				errmsg = msprintf(gettext("%s: Options for Gradient should be either ON or OFF"), "fminunc");
					error(errmsg);     	
        			end
        	case "Hessian" then
        			if (param(2*i)=="ON") then
        				//To check whether the user has provided Hessian function if Hessian Option is "ON"
						if (flag1==0) then
							if (rhs<4) then    
				     			errmsg = msprintf(gettext("%s: Hessian function is missing"), "fminunc");
				     			error(errmsg);     			
        					end
        					//This flag is activated(ie. =1) if Hessian is supplied
        					flag2=1;
        			        pos_h=4;
        					_h=varargin(4);
        				elseif (flag1==1) then
							if (rhs<5) then    
				     			errmsg = msprintf(gettext("%s: Hessian function is missing"), "fminunc");
				     			error(errmsg);     			
        					end
        					//This flag is activated(ie. =1) if Hessian is supplied
        					flag2=1;
        			        pos_h=5;
        					_h=varargin(5);
        				end
        			//To check whether Wrong entry(other than ON/OFF) is entered	            
        			elseif (param(2*i)~="ON" & param(2*i)~="OFF") then    
        				errmsg = msprintf(gettext("%s: Options for Hessian should be either ON or OFF"), "fminunc");
					error(errmsg);   
        			end
    		else
    	      		errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "fminunc", param(2*i-1));
    	      		error(errmsg)
    		end
   	end
   
   
	//Defining a function to calculate Gradient or Hessian if the respective user entry is OFF 
   	function y=_gradhess(x,t)
		if t==1 then	//To return Gradient
			y=numderivative(_f,x)		
		else		//To return Hessiam]n
			[grad,y]=numderivative(_f,x)
		end
   	endfunction
   	
   	//To check the correct no. of inputs given by the user
   	if (flag1==0 & flag2==0)
   		if(rhs>3) then
        	errmsg = msprintf(gettext("%s: Only 3 Inputs are Needed for this option(GradObj=OFF, HessObj=OFF), but %d were recorded"), "fminunc",rhs);
			error(errmsg); 
		end
   elseif ((flag1==1 & flag2==0) | (flag1==0 & flag2==1)) then
  		if(rhs>12) then
        	errmsg = msprintf(gettext("%s: Only 4 Inputs were needed for this option, but %d were recorded"), "fminunc",rhs);
			error(errmsg);
		end
   elseif (flag1==1 & flag2==1)
   		if(rhs>14) then
        	errmsg = msprintf(gettext("%s: Only 5 Inputs are Needed for this option(GradObj=ON, HessObj=ON), but %d were recorded"), "fminunc",rhs);
			error(errmsg); 
		end
	end
	
  //To check the correct input of Gradient and Hessian Functions from Users	     	
   if (flag1==1) then
   		if (type(_g) ~= 13 & type(_g) ~= 11) then
  			errmsg = msprintf(gettext("%s: Expected function for Gradient of Objective, since GradObj=ON"), "fminunc");
   			error(errmsg);
   		end
   		if(execstr('sample_g=_g(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Objective and x0 did not match "), "fminunc", rhs);
   			error(errmsg);
		end
		sample_g=_g(x0);
		if (size(sample_g,1)==s(2) & size(sample_g,2)==1) then
		elseif (size(sample_g,1)==1 & size(sample_g,2)==s(2)) then
		elseif (size(sample_g,1)~=1 & size(sample_g,2)~=1) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Gradient function(%dth Parameter)---->Row Vector function is Expected"), "fminunc",pos_g);
   			error(errmsg);
   		end
   	end
   	if (flag2==1) then
   		if (type(_h) ~= 13 & type(_h) ~= 11) then
  			errmsg = msprintf(gettext("%s: Expected function for Hessian of Objective, since HessObj=ON"), "fminunc");
   			error(errmsg);
   		end
   		if(execstr('sample_h=_h(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Hessian function of Objective and x0 did not match "), "fminunc", rhs);
   			error(errmsg);
		end
		sample_h=_h(x0);
   		if(size(sample_h,1)~=s(2) | size(sample_h,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Hessian function(%dth Parameter)---->Symmetric Matrix function is Expected "), "fminunc",pos_h);
   			error(errmsg);
   		end
   	end

    //Calling the Ipopt Function for solving the above Problem
	[xopt,fopt,status,iter,cpu,obj_eval,dual,gradient, hessian1] = solveminuncp(_f,_gradhess,flag1,_g,flag2,_h,x0,options);
   
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
    		gradient=[]
		hessian=[]
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
    if (flag1==1 |flag2==1) then
		disp("||||||Please Make sure you have entered Correct Functions for Gradient or Hessian -->Scilab Will Calculate Based on your input only||||||");
    end	
endfunction
