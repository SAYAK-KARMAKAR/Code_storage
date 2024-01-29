clc
clear
function[count,rmse,phi]=gseidel(A,r,phio,eps_max,omega)
    //Linear Equation:A*phi=r

[nr_A,nc_A]=size(A)//size of A
[nr_r,nc_r]=size (r)//size of r

if nc_A<>nr_A then
    error('A is not a square matrix')
    abort;
end
if nc_A<>nr_r then
    error('Not compatible matrices')
    abort;
end
n=nc_A

count=0
rmse=1
phi=phio
while rmse>eps_max
    rmse=0
    for i=1:1:n
        resi=r(i)
        for j=1:n
            resi=resi-A(i,j)*phi(j)
        end
        phi(i)=phi(i)+omega*resi/A(i,i)
        rmse=rmse+(omega*resi/A(i,i)).^2
    end
    rmse=sqrt(rmse/n)
    count=count+1;
    disp("Count rmse",[count rmse])
end
endfunction
//....Problem dependent parameters
nnode=31; //number of nodes
xl=0;//coordinate of left end boundary xL
xr=1000;
cconf=1e-11;
T=2e-5;
hs=90;
c0=90;
c1=0.06;
c2=-0.00003;
x=linspace(xl,xr,nnode);//In between left and right boundary, nnode numbers of x coordinates has been generated.
delx=x(2)-x(1);//gird size
//.....Initialization of matrices...
h=zeros(nnode,1);
A=zeros(nnode,nnode);
r=zeros(nnode,1);
//Left boundary
A(1,1)=1.0;
r(1)=hs;
//Interior nodes
for i=2:nnode-1
    A(i,i-1)=(1.0/(delx^2))*delx^2;
    A(i,i)=-((cconf/T)+2.0/(delx^2))*delx^2;
    A(i,i+1)=(1.0/(delx^2))*delx^2;
    r(i)=-(cconf/T)*(c0+c1*x(i)+c2*x(i)^2)*delx^2;//delx^2 is multiplied with each term for scaling
end
//Right boundary
//Two point
A(nnode,nnode)=(1/delx)*delx;
A(nnode,nnode-1)=-(1/delx)*delx;//delx is multiplied for scaling
r(nnode)=0;
//....
ho=hs*ones(nnode,1); //Initial guess =value specified at left hand boundary
eps_max=1e-6
omega=1
[count,rmse,h]=gseidel(A,r,ho,eps_max,omega)
plot(x,h','-g')
[L,U,E]=lu(A)
D=diag(diag(A))
for i=1:nnode
    L(i,i)=0
    U(i,i)=0
end
c=-inv(D)*(L+U)
eigv=max(abs(spec(c)))
disp("Eigenvalue",eigv)

