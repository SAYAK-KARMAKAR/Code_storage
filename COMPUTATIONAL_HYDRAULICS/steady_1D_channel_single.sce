clc
clear
//~~~~~~given data~~~~~~~~
Q=20; //m3/s
S0=0.0008;
n=0.015;
B=15;//m
g=9.81;//m/s^2
Lx=200 //m
yd=0.6;//m
mnode=201;
//no of nodes in x direction
eps_max=1e-6;
global('Q','S0','n','B','g')

//~~~~~~~~problem dependent parameters~~~~~~
alpha=1;
xc=linspace(0,Lx,mnode);
delta_x=Lx/(mnode-1);
zv(mnode)=0;
//elevation at the end section is considered to be at datum.

for i=mnode-1:-1:1
    zv(i)=zv(i+1)+S0*delta_x
    //downsream section was considered to be the end section which has lowest elevation. Elv. linealy increases in u/s direction.
end

yv=zeros(mnode,1);
 //It is the variable (water depth) defined for (Nl+1) number of nodes.
C1=alpha*Q^2/(2*g);
C2=(1/2)*n^2*Q^2*delta_x;
//c1 and c2 are constants. See equation (117)

global('C1','C2','delta_x')

function Av=areav(y)
    Av=B*y;
    //Cross-section area
endfunction

function dAv=dareav(y)
    dAv=B;
    //dAv gives dA/dy.
endfunction

function Rv=HRv(y)
    Rv=B*y/(B+2*y);
    //'Rv' gives hydraulic radius.
endfunction

function dRv=dHRv(y)
    dRv=B^2/(B+2*y)^2;
    //'dRv' gives dR/dy.
endfunction

//...........................
function Mliv=Mli(y1,y2)
    Mliv=(y2-y1)-S0*delta_x+C1*(areav(y2)^(-2)-areav(y1)^(-2))+C2*(HRv(y2)^(-4/3)*areav(y2)^(-2)+HRv(y1)^(-4/3)*areav(y1)^(-2));
    //It is momentum function M_{l,i}."-S0*delta_x" gives the the difference of the elevations of the two sections i.e."zv(2)-zv(1)". See equation (117).
endfunction

function dMdyiv=dMdyi(y)
    term1=(2/areav(y)^3*dareav(y));
    term2=2*areav(y)^(-3)*HRv(y)^(-4/3)*dareav(y);
    term3=(4/3)*areav(y)^(-2)*HRv(y)^(-7/3)*dHRv(y);
    dMdyiv=1+C1*term1-C2*(term2+term3);
    //dMdyiv is nothing but dM_{l,i}/dy_{l,i}. See equation (118)
endfunction

function dMdyip1v=dMdyip1(y)
    term1=(2/areav(y)^3)*dareav(y);
    term2=2*areav^(-3)*HRv^(-4/3)*dareav(y);
    term3=(4/3)*areav(y)^(-2)*HRv^(-7/3)*dHRv(y);
    dMdyip1v=1-C1*term1-C2*(term2+term3);
    //This gives dM_{l,i}/dy_{l,i+1}. See equation(119).
endfunction
A=zeros(mnode,mnode)
//Elements of jacobian matrix should be inserted in this 'A' matrix.
r=zeros(mnode,1);
count=0;
rmse=1;
yv=yd*ones(mnode,1);
//Downstream end value (yd) is taken as the initial guess value.

//~~~~~~~~~~~Space loop~~~~~~~~~~~~~
while rmse>eps_max 
    //While rmse>eps_max, then only we should iterate. 
rmse=0;
for i=1:mnode-1
    A(i,i)=dMdyi(yv(i));
    A(i,i+1)=dMdyip1(yv(i+1)); 
    r=Mli(yv(i),yv(i+1));
end
//~~~~~Subcritical boundary conditions~~~~~
A(mnode,mnode)=1;
r(mnode)=-(yv(mnode)-yd);
//Here, 'r' vector in "[A]{dely}={r}" equation contains "negative of momentum function"(see equation 123). From 116.1, we have d/s boundary condition i.e., "DB_{l,Nl+1}=y_{l,Nl+1}-y_{d}". Hence, value of 'r(mnode)' is justified.

dely=A\r;
for i=1:mnode
    yv(i)=yv(i)+dely(i);
    // Initially, yv(i)=yd; for all values of i. Then, it got modified. 
    rmse=rmse+dely(i)^2;
    //Actually, using this loop, we are summing up the square of all the errors. RMSE was calculated later.( I.T--> initial rmse should be considered as 0. why 1?) 
end

rmse=sqrt(rmse/mnode);
count=count+1;
//disp('COUNT RMSE')
//disp([count; rmse])
end

//~~~~~~figures and plots~~~~~~~~
 plot(xc',yv+zv,"-r")
 plot([0 200],[zv(1) zv(mnode)],'b-')
 xtitle("steady single channel flow","x axis","Flow depth")
