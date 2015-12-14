function y=f(x)
y=3*x(1)^2+2*x(1)*x(2)+x(2)^2-4*x(1)+5*x(2)
endfunction
pt=[1,1]
exec builder.sce
exec loader.sce
[x,f,e,s,g,h]=fminunc(f,pt)

