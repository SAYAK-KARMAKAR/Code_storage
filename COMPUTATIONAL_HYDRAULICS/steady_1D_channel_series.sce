clc
clear
//Question is at(5 no.-p-24). Lecture 37
//This problem considers constant discharge Q. flow depth 'y' is the only variable here.
//~~~~~~Given Data~~~~~~
chl=2;
//We have two channels in series.
Q=20;//m^3/S
S0=[0.0004 0.0008];
n=[0.010 0.015]
B=15;//m
g=9.81;//m\s^2
Lx=[100 100];
yd=0.6;//m
mnode=[101 101];
eps_max=1e-6;
global('Q','B','g')
//Here, S0 and n are varying. So, we don't keep those as global variables.

//~~~~~problem dependent parameters~~~~~~~
alpha=[1,1];
yv=zeros(sum(mnode),1);
// Here,sum(mnode)=sum of all elements in 'mnode' matrix.Therefore, yv=zeros(202,1). Because, total number of sections=101*2=202.
for l=1:chl
//Here, l is the channel numbering. chl is number of channels(=2 here).
    delta_x(l)=Lx(l)/(mnode(l)-1);
    C1(l)=alpha(l)*Q^2/(2*g);
    C2=(1/2)*n(l)^2*Q^2*delta_x(l);
    //Here, two loops will run.  for l=1 and l=2.
end
mc=sum(mnode);
// It gives total number of sections (i.e.sum of all elements in 'mnode' matrix). Because, we need to get those many 'depth' values.


for l=chl:-1:1
// i.e. l=2:-1:1 ,for this case.

    for i=mnode(l):-1:1
        if(l==chl & i==mnode(l))then
            zv(mc)=0;
 //Here, 'l' represents channel numbering. Number of channel is 2. So, chl=2. Therefore, this condition means l==2 and i=mnode(2)=101; i.e. the downstream section.
        end
        if(l<>chl & i==mnode(l)) then
// end section of a channel other than the last channel.
            mc=mc-1;
            zv(mc)=zv(mc+1)
// This 'if' loop will be executed (i.e. will take a value ,other than zero, for 101-th row of the 'zv' vector) when the next 'if' loop will evaluate the zv value of the 1st section of the second channel(i.e. 102-th entry of 'zv' vector).
        end
        if (i<>mnode(l)) then
            mc=mc-1; 
            zv(mc)=zv(mc+1)+S0(l)*delta_x(l);
//This 'if' loop is for all other general setions. In first iteration, this 3rd 'if' loop will evaluate nothing. Only 1st 'if' loop will be executed for l=2 & i==mnode(2)=101-->zv=0. for 2nd loop, 
//l=2,i=mnode-1=100 and mc=201; (i.e. zv(201)=zv(202)+S0(2)*delta_x(2).
//Then, l=2,i=99 and mc=200; (i.e. zv(200)=zv(201)+S0(2)*delta_x(2)...and so on.....
//When "l=2, i=1 and mc=102; (i.e. zv(102)=zv(103)+S0(2)*delta_x(2)" is done, then 2nd if loop will take non-zero value.i.e. l=1,i=101 & mc=102-1=101--> zv(101)=zv(101+1)=zv(102). Because last section of a channel and the 1st section of the very next channel will have same elevation. 
        end
    end  
end
xv=[linspace(0,Lx(1),mnode(1)) linspace(Lx(1),Lx(1)+Lx(2),mnode(2))]

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
function Mliv=Mli(y1,y2,S0,delta_x,C1,C2)
//Previously, "S0,delta_x,C1,C2" were global values. But, now these are different for different channels.
    Mliv=(y2-y1)-S0*delta_x+C1*(areav(y2)^(-2)-areav(y1)^(-2))+C2*(HRv(y2)^(-4/3)*areav(y2)^(-2)+HRv(y1)^(-4/3)*areav(y1)^(-2));
endfunction

function dMdyiv=dMdyi(y,C1,C2)
    term1=(2/areav(y)^3*dareav(y));
    term2=2*areav(y)^(-3)*HRv(y)^(-4/3)*dareav(y);
    term3=(4/3)*areav(y)^(-2)*HRv(y)^(-7/3)*dHRv(y);
    dMdyiv=1+C1*term1-C2*(term2+term3);
    //dMdyiv is nothing but dM_{l,i}/dy_{l,i}. See equation (118)
endfunction

function dMdyip1v=dMdyip1(y,C1,C2)
    term1=(2/areav(y)^3)*dareav(y);
    term2=2*areav(y)^(-3)*HRv^(-4/3)*dareav(y);
    term3=(4/3)*areav(y)^(-2)*HRv^(-7/3)*dHRv(y);
    dMdyip1v=1-C1*term1-C2*(term2+term3);
    //This gives dM_{l,i}/dy_{l,i+1}. See equation(119).
endfunction

A=zeros(sum(mnode),sum(mnode))
//Elements of jacobian matrix should be inserted in this 'A' matrix.No of rows and columns in A matrix= total number of nodes/sections in all channels.
r=zeros(sum(mnode),1);
count=0;
rmse=1;
yv=yd*ones(sum(mnode),1);
//Downstream end value (yd) is taken as the initial guess value.

//~~~~~~~~~~~Space loop~~~~~~~~~~~~~
while rmse>eps_max 
rmse=0;
mc=0
//Here, mc is nothing but the counter of iteration. 
for l=1:chl
    for i=1:mnode(l) 
        mc=mc+1
        if (l==chl & i==mnode(l)) then
            A(mc,mc)=1;
            r(mc)=-(yv(mc)-yd);
//This if loop is applicable for downstream section of the last channel. For this 'if' loop,  mc=0+1=1. if l=2 & i=mnode(2)=101 then, A(202,202)=1; r(202)=-(yv(1)-yd).
// This satisfy the equation "DB_{2,mnode}=yv_{l,Nl+1}-y_{d}"(see equation 116.1).So, r(202)=-DB_{2,mnode}=-(yv(1)-yd);(See equation 123).
//Clearly, this 'if' loop will be executed at last iteration for "l=chl & i=mnode(= 101)".
        else
            if (i==mnode(l)) then
                A(mc,mc)=1;
                A(mc,mc+1)=-1;
                r(mc)=-(yv(mc)-yv(mc+1));
//This 'if' loop is applicable for the last section of all channel(s) other than end one.For this problem, this 'if' loop will be executed for 101-th iteration i.e  when mc=101.
//Here, for l=1, if i==mnode(1)=101,mc=101--> A(101,101)=1, A(101,102)=-1, r(101)=-(yv(101)-yv(102))
//This should satisfy the junction condition (Depth continuity) "y_{1,mnode}=y_{2,1}". Therefore in equation form,"Mliv=y_{1,mnode}-y_{2,1}". So,r(101)=-Mliv=-(y_{1,mnode}-y_{2,1}).(See equation 123).
        else
            A(mc,mc)=dMdyi(yv(mc),C1(l),C2(l));
            A(mc,mc+1)=dMdyip1(yv(mc+1),C1(l),C2(l));
            r(mc)=-Mli(yv(mc),yv(mc+1),S0(l),delta_x(l),C1(l),C2(l));
        end
    end
    end             
end

dely=A\r;
for i=1:sum(mnode)
    yv(i)=yv(i)*dely(i);
    rmse=rmse+dely(i)^2;
end
rmse=sqrt(rmse/sum(mnode));
count=count+1;
disp([count rmse])
end
//Figures and Plots
plot(xv,yv+zv,"-r")
plot([Lx(1) Lx(1)+Lx(2)],[zv(mnode(1)+1) zv(sum(mnode))],'m-')
plot([0,Lx(1)],[zv(1),zv(mnode(1))],'m-')
xtitle('steady channel flow', 'x axis','Flow depth')
//zv(1)
//zv(mnode(1))




































