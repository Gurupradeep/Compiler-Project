function y=f(x)
y=x(1)^2+x(2)^2
endfunction
pt=[1,2]
exec builder.sce
exec loader.sce
[x,f,e,s,g,h]=fminunc(f,pt)

