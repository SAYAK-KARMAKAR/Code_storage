clc
clear
function phi=ludcomp(A,r)
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
//Decomposition
for k=1:1:n-1
    for i=k+1:1:n
        gam=A(i,k)/A(k,k)
        A(i,k)=gam
        for j=k+1:n
            A(i,j)=A(i,j)-gam*A(k,j)
        end
    end
end

//Forward substitution
psi(1)=r(1)
for i=2:1:n
    sumj=r(i)
    for j=1:i-1
        sumj=sumj-A(i,j)*psi(j)
    end
    psi(i)=sumj
end

//Backward substitution
phi(n)=psi(n)/A(n,n)
for i=n-1:-1:1
    sumj=psi(i)
    for j=i+1:n
        sumj=sumj-A(i,j)*phi(j)
    end
    phi(i)=sumj/A(i,i)
end
endfunction
//.....
function FV=Fun(phi)
    FV(1)=phi(1)^2+phi(1)*phi(2)-3
    FV(2)=phi(1)+phi(2)^3-phi(2)*phi(3)-3
    FV(3)=-phi(2)^2+phi(3)^2+phi(3)*phi(4)-17
    FV(4)=-phi(3)^2+phi(4)^2+phi(4)*phi(5)-27
    FV(5)=phi(4)+phi(5)^2-29
endfunction
function JV=JM(phi)
    JV=zeros(length(phi),length(phi))
    //1
    JV(1,1)=2*phi(1)+phi(2)
    JV(1,2)=phi(1)
    //2
    JV(2,1)=1
    JV(2,2)=3*phi(2)^2-phi(3)
    JV(2,3)=-phi(3)
    //3
    JV(3,2)=-2*phi(2)
    JV(3,3)=2*phi(3)+phi(4)
    JV(3,4)=phi(3)
    //4
    JV(4,3)=-2*phi(3)
    JV(4,4)=2*phi(4)+phi(5)
    JV(4,5)=phi(4)
    //5
    JV(5,4)=1
    JV(5,5)=2*phi(5)
endfunction
function[count,rmse,phi]=newtonrn(phio,n,eps_max)
    count=0
    rmse=1
    while rmse>eps_max
        rmse=0
        JV=JM(phio)
        FV=Fun(phio)
        dphi=ludcomp(JV,-FV)
        for i=1:n
        phi(i)=phio(i)+dphi(i)
        rmse=rmse+(dphi(i)).^2
    end
    phio=phi
    rmse=sqrt(rmse/n)
    count=count+1
end
endfunction
              
//
phio1=[1
1
1
1
1];
phio2=[0
0
0
0
0];
phio3=[10
1
100
6
11];
eps_max=1e-6;
n=5
[count,rmse,phi]=newtonrn(phio3,n,eps_max)
FV=Fun(phi)
disp(phi)
    
    
    
    
