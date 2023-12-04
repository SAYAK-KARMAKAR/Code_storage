//GRADUALLY VARIED FLOW- EULER CAUCHY METHOD
clc
clear all
//Data given
Q=20;
S0=0.0008;
B=15;//m, Channel width @rectangular channel
Lx=200;//length of the channel reach
y0=0.8;//initial depth at x=0
n=0.015;
g=9.81;
mnode=201;//number of nodes in the x directions
global('Q','S0','n','B','g')//global variables . for writing multiple functions later, we can directly utilize these constant values.

//***********************************************
//calculation of critical depth
yc=(Q^2/(g*B^2))^(1/3);//for rectangular channel.equation (73)
disp("Critical depth",yc)
//***************************************************
//calculation of normal depth using NR method [from equation (75.1)]

//To define G and G' functions
function Gval=Gfunc(y)
    Gval=(S0^(1/2)*B^(5/3)/n)*(y/(B+2*y))^(2/3)*y-Q;//eqn(74)
endfunction
function Gpval=Gderi(y)
    term1=(S0^(1/2)*B^(5/3))/(3*n);
    term2=y^(2/3)*(5*B+6*y);
    term3=(B+2*y)^(5/3);
//    Gpval=(S0^(1/2)*B^(5/3))/(3*n)*y^(2/3)*(5*B+6*y)/(B+2*y)^(5/3);//G' from eqn (76)
Gpval=term1*term2/term3;
endfunction

//Newton-Raphson method (to solve nonlinear 'Gval' eqaution for getting yn)
eps_max=1e-6;
aerror=1; //absolute error
yn=yc; //To start the iteration loop, critical depth is our initial guess.
while abs(aerror)>eps_max
    aerror=(Gfunc(yn)/Gderi(yn));
    yn=yn-aerror; //yn(p)=yn(p-1)-(G(yn(p-1))/G'(yn(p-1)))..eqn(75.1))
end
disp("Normal Depth",yn)
//****************************************

//dydx function calculation from Governing Equation
function dydx=psi(x,y)  //equation (73)
    A_y=B*y;
    P_y=B+2*y; 
    R_y=A_y/P_y;//Area, Perimeter and hydraulic radius are function of y
    Sf=(n^2*Q^2)/(R_y^(4/3)*A_y^2);//friction slope
    Frs=(Q^2*B)/(g*A_y^3);//Fr^2
    dydx=(S0-Sf)/(1-Frs);
endfunction
//*********************************************
//Euler-Cauchy method [to get GVF profile for nodes along x direction]
xc=linspace(0,Lx,mnode);//x coordinate vector of nodes
delta_x=Lx/(mnode-1);
yv=zeros(mnode,1);//Initialization of flow depth vector
yv(1)=y0;//Input Boundary value in yv vector
//input other entries
for i=2:mnode
    K1=delta_x*psi(xc(i-1),yv(i-1)); //K1 and K2 for ith point are calculated from previous (i-1)th point data.(explicit) 
    K2=delta_x*psi(xc(i-1)+delta_x,yv(i-1)+K1)
    
    yv(i)=yv(i-1)+0.5*(K1+K2);//Equation (78)
end

plot(xc,yv,"-b")
set(gca(),"auto_clear","off")
plot([0 Lx],[yn yn],'-m')
set(gca(),"auto_clear","off")
plot([0 Lx],[yc yc],'-b')
xtitle("GVF Euler Cauchy Method","x axis(m)","flow depth")







 
