//GRADUALLY VARIED FLOW-IMPLICIT RK2
clc
clear 
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
disp("yc",yc)
//***************************************************
//calculation of normal depth using NR method [from equation (75.1)]

//To define G and G' functions
function Gval=Gfunc(y)
    Gval=(S0^(1/2)*B^(5/3)/n)*(y/(B+2*y))^(2/3)*y-Q;//eqn(127)
endfunction
function Gpval=Gderi(y)
    term1=(S0^(1/2)*B^(5/3))/(3*n);
    term2=y^(2/3)*(5*B+6*y);
    term3=(B+2*y)^(5/3);
//    Gpval=(S0^(1/2)*B^(5/3))/(3*n)*y^(2/3)*(5*B+6*y)/(B+2*y)^(5/3);//G' 
Gpval=term1*term2/term3;//from eqn (129)
endfunction

//Newton-Raphson method (to solve nonlinear 'Gval' eqaution for getting yn)
eps_max=1e-6;
aerror=1; //absolute error
yn=yc; //To start the iteration loop, critical depth is our initial guess.
while abs(aerror)>eps_max
    aerror=(Gfunc(yn)/Gderi(yn));
    yn=yn-aerror; //yn(p)=yn(p-1)-(G(yn(p-1))/G'(yn(p-1)))..eqn(75.1))
end
disp("yn",yn)
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
//function Fval=Ffunc(y,x,delta_x,yp)
//    Fval=y-delta_x*psi(x,y)-yp //y=new value,yp=previous known value.Equation(134)
//endfunction
//function Fpval=Fderi(y,delta_x) //derivative of Fval wrt y(new)
//    term1=(1-Q^2/(B^2*g*y^3))^(-1);
//    term2=(2*n^2*Q^2)/(B^2*y^3);
//    term3=(B*y/(B+2*y))^(-4/3);
//    term4=(4*n^2*Q^2)/(3*B^2*y^2);
//    term5=(B*y/(B+2*y))^(-7/3);
//    term6=((B/(B+2*y))-(2*B*y/(B+2*y)^2));
//    term7=(3*Q^2)/(B^2*g*y^4);
//    term8=(1-Q^2/(B^2*g*y^3))^(-2);
//    term9=S0;
//    term10=(n^2*Q^2)/(B^2*y^2);
//    term11=(B*y/(B+2*y))^(4/3);
//    Fpval=1-delta_x*(term1*(term2*term3+term4*term5*term6)-term7*term8*(term9-term10*term11));   
//endfunction
/////function psipval=psideri(x,y,delta_x) 
function psipval=psideri(x,y)//derivative of psival wrt y(new)
    term1=(1-Q^2/(B^2*g*y^3))^(-1);
    term2=(2*n^2*Q^2)/(B^2*y^3);
    term3=(B*y/(B+2*y))^(-4/3);
    term4=(4*n^2*Q^2)/(3*B^2*y^2);
    term5=(B*y/(B+2*y))^(-7/3);
    term6=((B/(B+2*y))-(2*B*y/(B+2*y)^2));
    term7=(3*Q^2)/(B^2*g*y^4);
    term8=(1-Q^2/(B^2*g*y^3))^(-2);
    term9=S0;
    term10=(n^2*Q^2)/(B^2*y^2);
    term11=(B*y/(B+2*y))^(4/3);
    psipval=(term1*(term2*term3+term4*term5*term6)-term7*term8*(term9-term10*term11));//using eqn(135.1)
    endfunction
//*********************************************
//Implicit Backward Euler method [to get GVF profile for nodes along x direction]
xc=linspace(0,Lx,mnode);//x coordinate vector of nodes
delta_x=Lx/(mnode-1);
yv=zeros(mnode,1);//Initialization of flow depth vector
yv(1)=y0;//Input Boundary value in yv vector
for i=2:mnode
    K1=delta_x*(1-0.5*delta_x*psideri(xc(i-1)+0.5*delta_x,yv(i-1)))^(-1)*psi(xc(i-1)+0.5*delta_x,yv(i-1));//using semi-implicit eqn (138) and putting values from (135.1) and (136.1)
    yv(i)=yv(i-1)+1*K1;
end
//plot
plot(xc,yv,"-g")
set(gca(),"auto_clear","off")
plot([0 Lx],[yn yn],'-m')//Normal depth line
set(gca(),"auto_clear","off")
plot([0 Lx],[yc yc],'-b') //Critical depth line
xtitle("GVF implicit RK2","x axis(m)","flow depth(m)")







 
