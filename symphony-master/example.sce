function y=fun(x)
y=3*x(1)^2+2*x(1)*x(2)+x(2)^2-4*x(1)+5*x(2)
endfunction
pt=[1,1]
exec builder.sce
exec loader.sce
options=list("MaxIter", [4],"CpuTime", [100])
[x,f,e,s,g,h]=fminunc(fun,pt,options)

