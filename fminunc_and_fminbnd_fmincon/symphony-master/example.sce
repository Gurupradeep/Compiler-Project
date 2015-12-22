//fminunc
function y=fun(x)
y=100*(x(2)-x(1)^2)^2 + (1-x(1))^2
endfunction

function y=grad(x)
//y=[3*x(1)^2,3*x(2)^2];
y= [-400*x(1)*x(2) + 400*x(1)^3 + 2*x(1)-2, 200*(x(2)-x(1)^2)];
endfunction

function y=hess(x)
//y=[6*x(1),0;0,6*x(2)]
y= [1200*x(1)^2- 400*x(2) + 2, -400*x(1);-400*x(1), 200 ];
endfunction

//FAILURE CASES
//Fails sometimes if starting point is a stationary point i.e. f'=0. (works for x1^2+x2^2 and (0,0) but fails for (x1-1)^2+(x2-1)^2 and (1,1))
//Fails when it converges to point of inflecion. So if function has a point of inflection nearby when compared to the local minimum the  function may fail to find the optimal value.

pt=[1,1.1];
options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", "OFF", "Hessian", "OFF");
[x,f,e,s,g,h]=fminunc(fun,pt,options)


//fminbnd
function y=fun1(x)
y=(x-1)^3
endfunction

//FAILURE CASES
//Fails if point of inflection is x=0 and is included in the bounded interval

options1=list("MaxIter",[15000000],"CpuTime", [10000])
[x1,f1,e1,s1]=fminbnd(fun1,-10,12,options1)

