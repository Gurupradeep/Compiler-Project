// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: harpreet.mertia@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


//function [xopt,fopt,exitflag,output,lambda] = qpipopt (varargin)
  function [xopt,fopt,exitflag,output,gradient,hessian] = fminunc (varargin)

//To check the number of input and output argument
   [lhs , rhs] = argn();
	
//To check the number of argument given by user
   if ( rhs<2 | rhs>3 ) then
    errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 2 or 3"), "fminunc", rhs);
    error(errmsg)
   end
   
   
   fun_ = varargin(1);
   x0 = varargin(2);
   y=fun_(x0);
   
   if ( rhs<3 | size(varargin(3)) ==0 ) then
      param = list(); 
   else
      param =varargin(3);
   end
   

   options = list("MaxIter", [3000],"CpuTime", [600]);
      

   for i = 1:(size(param))/2
       	select param(2*i-1)
    	case "MaxIter" then
          		options(2*i) = param(2*i);
       	case "CpuTime" then
          		options(2*i) = param(2*i);
    	else
    	      errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "fminunc", param(2*i-1));
    	      error(errmsg)
    	end
   end


	function y=gradhess_(x,t)
		if t==1 then
			y=numderivative(fun_,x)
		else
			[grad,y]=numderivative(fun_,x)
		end
	endfunction

   [xopt,fopt,status,iter,gradient, hessian1] = solveminuncp(fun_,gradhess_,x0,options);
   xopt = xopt';
   exitflag = status;
   output = struct("Iterations"      , []);
   output.Iterations = iter;

    s=size(gradient)
    for i =1:s(2)
    	for j =1:s(2)
		hessian(i,j)= hessian1(j+((i-1)*s(2)))
	end
    end

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
    

endfunction
