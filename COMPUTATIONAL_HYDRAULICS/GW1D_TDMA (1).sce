//Groundwater 1D using Tridiagonal matrix method
clc
clear
function phi=tdmasolv(b,d,a,r)
    
//n:number of rows
n=length(d);

a(1)=a(1)/d(1);
r(1)=r(1)/d(1);
//forward elimination
for i=2:n-1
    fact=d(i)-b(i)*a(i-1);
    a(i)=a(i)/fact;
    r(i)=(r(i)-b(i)*r(i-1))/fact;
end

r(n)=(r(n)-b(n)*r(n-1))/(d(n)-b(n)*a(n-1));

//Backward substitution
phi(n)=r(n);
for i=n-1:-1:1
    phi(i)=r(i)-a(i)*phi(i+1);
end
endfunction
//....Problem dependent parameters
nnode=21; //number of nodes
xl=0;//coordinate of left end boundary xL
xr=1000;
cconf=1e-11;
T=2e-5;
hs=90;
c0=90;
c1=0.06;
c2=-0.00003;
//....
x=linspace(xl,xr,nnode);//In between left and right boundary, nnode numbers of x coordinates has been generated.
delx=x(2)-x(1);//gird size
//.....Initialization of matrices...
h=zeros(nnode,1);
a=zeros(nnode,1);
d=zeros(nnode,1);
b=zeros(nnode,1);
r=zeros(nnode,1);
//Left boundary
d(1,1)=1.0;
r(1)=hs;
//Interior nodes
for i=2:nnode-1
    b(i,1)=1.0/(delx^2);
    d(i,1)=-((cconf/T)+2.0/(delx^2));
    a(i,1)=1.0/(delx^2);
    r(i)=-(cconf/T)*(c0+c1*x(i)+c2*x(i)^2);
end
//Right boundary
//Two point
d(nnode,1)=1/delx;
b(nnode,1)=-1/delx;
r(nnode)=0;

h=tdmasolv(b,d,a,r)
plot(x,h','-r')
disp("Head at nodes",h)
disp("maximum head",max(h))
