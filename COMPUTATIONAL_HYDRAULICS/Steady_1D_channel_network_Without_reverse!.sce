clc
clear
function Av=areav(y,B,m1,m2)
    Av=B*y+(1/2)*(m1+m2)*y^2;
endfunction
function dAv=dareav(y,B,m1,m2)
    dAv=B+(m1+m2)*y;
//dAv means dA/dy. See (157.1).
endfunction
function Rv=HRv(y,B,m1,m2)
    Rv=(B*y+(1/2)*(m1+m2)*y^2)/(B+(sqrt(1+m1^2)+sqrt(1+m2^2))*y);
// it is hydraulic radius. see equation (155)
endfunction
function dRv=dHRv(y,B,m1,m2)
    Tw=B+(m1+m2)*y;
//Top width. See equation (155)
    Pm=(B+(sqrt(1+m1^2)+sqrt(1+m2^2))*y);
//Wetted perimeter. See equation (155)
    Rh=HRv(y,B,m1,m2);
    dPdy=(sqrt(1+m1^2)+sqrt(1+m2^2));
//dP/dy. See equation (157.2).
    dRv=(Tw/Pm)-(Rh/Pm)*dPdy;
//dRv means dR/dy. See equation (120)for its derrivation.    
endfunction
function Mliv=Mli(y1,Q1,y2,Q2,S0,delta_x,D1,D2,B,m1,m2)
    Mliv=(y2-y1)-S0*delta_x+D1*(Q2^2*areav(y2,B,m1,m2)^(-2)-Q1^2*areav(y1,B,m1,m2)^(-2))+D2*(Q2*abs(Q2)*HRv(y2,B,m1,m2)^(-4/3)*areav(y2,B,m1,m2)^(-2)+Q1*abs(Q1)*HRv(y1,B,m1,m2)^(-4/3)*areav(y1,B,m1,m2)^(-2)); 
endfunction

function dMdyiv=dMdyi(y,Q,D1,D2,B,m1,m2)
    term1=(2*Q^2/areav(y,B,m1,m2)^(3))*dareav(y,B,m1,m2);
    term2=2*Q*abs(Q)*areav(y,B,m1,m2)^(-3)*HRv(y,B,m1,m2)^(-4/3)*dareav(y,B,m1,m2);
    term3=(4/3)*Q*abs(Q)*areav(y,B,m1,m2)^(-2)*HRv(y,B,m1,m2)^(-7/3)*dHRv(y,B,m1,m2);
    dMdyiv=-1+D1*term1-D2*(term2+term3);
    //From equation (147.1).Evaluates dM_{l,i}/dy_{l,i}.
endfunction

function dMdyip1v=dMdyi(y,Q,D1,D2,B,m1,m2)
    term1=(2*Q^2/area(y,B,m1,m2)^(3))*dareav(y,B,m1,m2);
    term2=2*Q*abs(Q)*areav(y,B,m1,m2)^(-3)*HRv(y,B,m1,m2)^(-4/3)*dareav(y,B,m1,m2);
    term3=(4/3)*Q*abs(Q)*areav(y,B,m1,m2)^(-2)*HRv(y,B,m1,m2)^(-7/3)*dHRv(y,B,m1,m2);
    dMdyiv=-1-D1*term1-D2*(term2+term3);
    //From equation (147.3). Evaluates dM_{l,i}/dy_{l,i+1}.
endfunction

function dMdQip1v=dMdQip1(y,Q,D1,D2,B,m1,m2)
    term1=2*Q*areav(y,B,m1,m2)^(-3);
    term2=2*abs(Q)*areav(y,B,m1,m2)^(-2)*HRv(y,B,m1,m2)^(-4/3);
    dMdQip1v=D1*term1+D2*term2;
    //From equation (147.4).
endfunction

function dMdQiv=dMdQi(y,Q,D1,D2,B,m1,m2)
    term1=2*Q*areav(y,B,m1,m2)^(-3);
    term2=2*abs(Q)*areav(y,B,m1,m2)^(-2)*HRv(y,B,m1,m2)^(-4/3);
    dMdQip1v=-D1*term1+D2*term2;
    //From equation (147.2).
endfunction

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Channel reach: start + , end - .
//Flow depth condition: 1
//Flow rate (discharge) condition: 2
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~Given Data~~~~~~~~~
g=9.81; //m/s^2
global('g')
yd=3; //m
Qd=250; //m^3/s
eps_max=1e-6;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
junn=3;//number of junctions
chln=4;//number of channels
//~~~~Chl | Length | Width | m1 | m2 | Segment | n | S0 | JN1 | JN2 |
chl_inf=[1  100 50 2 2 25 0.0120 0.0005 0 3
         2 1500 30 2 2 75 0.0125 0.0004 3 2
         3  500 20 2 2 25 0.0130 0.0012 3 2 
         4  100 40 2 2 25 0.0135 0.0005 2 1];
jun_inf=[yd -Qd
         -99999 -99999
         -99999 -99999];
jun_con= [1 -4 0 0
          3 4 -3 -2
          3 -1 2 3];
//In 'jun_con' matrix, 1st column represents number of channels connected to that junction.
//Positive sign means 1st section of that channel is connected to that junction. Negative sign means (Nl+1)th section of that channel is connected to that junction.    
alpha=[1 1 1 1];

//~~~~~~~Derived informations from 'chl_inf' matrix~~~~~~~~~~
Lx=chl_inf(1:chln,2);
//We know, 'chln' means channel number=4.
B=chl_inf(1:chln,3);
m1=chl_inf(1:chln,4);
m2=chl_inf(1:chln,5);
delta_x=chl_inf(1:chln,6);
n=chl_inf(1:chln,7);
S0=chl_inf(1:chln,8);

mnode=Lx./delta_x+1;
//Here, mnode is calculated at point by point basis. Here, Lx=[100; 1500; 500; 100], delta_x=[25;75; 25; 25] results mnode=[5; 21; 21; 5].

//~~~~~~~Problem Dependent Parameters~~~~~~~
yv=yd*ones(sum(mnode),1);
Qv=Qd*ones(sum(mnode),1);
//Problem is initialized with downstream depth and discharge values.
gv=zeros(2*sum(mnode),1);
//'gv' will store the general variable with y and Q.

//~~~General identification matrix~~~
idv=0;
for l=1:chln
    for i=1:mnode(l)
        idv=idv+1;
        gid(l,i)=idv;
//Number of rows in 'gid' matrix= no. of channels=4 and number of columns= maximum no. of sections in a channel reach=max(5,21,21,5)=21. So, 'gid' is a 4*21 matrix.
    end
end
for l=1:chln
     for i=1:mnode(l)
     gv(2*gid(l,i)-1)=yv(gid(l,i))
     gv(2*gid(l,i))=Qv(gid(l,i))
//Initial Values of 'yv' and 'Qv' matries are stored in 'gv' or general variable matrix togetherly.
     end
end

for l=1:chln
    D1(l)=alpha(l)/(2*g);
    D2(l)=(1/2)*n(l)^2*delta_x(l);
end

A=zeros(2*sum(mnode),2*sum(mnode));
//'A' will store coff(s) of jacobian matrix.
r=zeros(2*sum(mnode),1);

Count=0;
rmse=1;
//~~~~~~~~~~~~Space loop~~~~~~~~~~~~~~~~~~~
while rmse>eps_max
    rmse=0;
    eqn=0;//Equation number
    
    //~~~~~Equations corresponding to segments(2N1+2N2+2N3+2N4)~~~~~~
    for l=1:chln
    for i=1:mnode(l)-1        
//~~~~~~Jacobians for Continuity eqn~~~~~~~
        eqn=eqn+1;
        A(eqn,2*gid(l,i)-1)=0;
//del(C_{l,i})/del(y_{l,i}).See equation (156).
        A(eqn,2*gid(l,i))=-1;
//del(C_{l,i})/del(Q_{l,i}).
        A(eqn,2*gid(l,i+1)-1)=0;
//del(C_{l,i})/del(y_{l,i+1}).
        A(eqn,2*gid(l,i+1))=1;
//del(C_{l,i})/del(Q_{l,i+1}).
        r(eqn)=0;

//~~~~~~~Jacobians for Momentum eqn~~~~~~~~
        eqn=eqn+1;
        A(eqn,2*gid(l,i)-1)=dMdyi(yv(gid(l,i)),Qv(gid(l,i)),D1(l),D2(l),B(l),m1(l),m2(l));
//del(M_{l,i})/del(y_{l,i}).See page 37 & equation (147.1).
         A(eqn,2*gid(l,i))=dMdQi(yv(gid(l,i)),Qv(gid(l,i)),D1(l),D2(l),B(l),m1(l),m2(l));
//del(M_{l,i})/del(Q_{l,i}).See equation (147.2).
          A(eqn,2*gid(l,i+1)-1)=dMdyip1(yv(gid(l,i+1)),Qv(gid(l,i+1)),D1(l),D2(l),B(l),m1(l),m2(l));
//del(M_{l,i})/del(y_{l,i+1}).See equation (147.3).
          A(eqn,2*gid(l,i+1))=dMdQip1(yv(gid(l,i+1)),Qv(gid(l,i+1)),D1(l),D2(l),B(l),m1(l),m2(l));
//del(M_{l,i})/del(Q_{l,i+1}).See equation (147.4).
          
          r(eqn)=-Mli(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),S0(l),delta_x(l),D1(l),D2(l),B(l),m1(l),m2(l));
          end
end

//~~~~~~~Junction Continuity Conditions~~~~~~~
for j=1:junn
    eqn=eqn+1;
    if(jun_inf(j,2)<>-99999) then
        r(eqn)=-jun_inf(j,2);
//////How??????????? I.T--> r=-DB=-(Q_{4,N4+1}+Q_d)...See Equation(159). But here, r=-(-Q_d)=Q_d....
    else
        r(eqn)=0;
//'r(eqn)=0' -->it means, no quantity is added there.
end

//    jun_con= [1 -4 0 0
//              3 4 -3 -2
//              3 -1 2 3];
for l=1:jun_con(j,1)
//1st column gives the information about number of channel connected to particular junction.
if (abs(jun_con(j,l+1))>eps_max) then
//Checks whether a non-zero value is there or not?!zero value means that, no channel is connected.
    if (jun_con(j,l+1)>0) then
//Checks whether the value is positive or negative. If positive, then the 1st section of the channel is connected to that junction. If negative, then the end section of the channel is connected to that junction.
        jn_node=1;
        A(eqn,2*gid(abs(jun_con(j,l+1)),jn_node))=-1;
        r(eqn)=r(eqn)-Qv(gid(abs(jun_con(j,l+1)),jn_node));
    end
    if (jun_con(j,l+1)<0) then
            jn_node=mnode(abs(jun_con(j,l+1)));
            A(eqn,2*gid(abs(jun_con(j,l+1)),jn_node))=1;
            r(eqn)=r(eqn)+Qv(gid(abs(jun_con(j,l+1)),jn_node));
    end
end
end
   r(eqn)=-r(eqn); 
end
//Explanation of the junction continuity condition: 
/////For 1st iteration: for j=1;
// eqn=eqn+1; 
//jun_inf(j,2)=jun_inf(1,2)=-250<>-99999;
// r(eqn)=-jun_inf(j,2)=-jun_inf(1,2)=-(-250)=250;
// l=1(i.e. for j=1 & l=1),
//if abs(jun_con(j,l+1))=abs(jun_con(1,2))=abs(-4)=4 > eps_max;
// if (jun_con(j,l+1)=jun_con(1,2)=-4<0)
// then jn_node=mnode(abs(jun_con(j,l+1)))=mnode(abs(jun_con(1,2)))=mnode(4)=5;
//A(eqn,2*gid(abs(jun_con(j,l+1)),jn_node))=A(eqn,2*gid(4,5))=A(eqn,2*52)=A(eqn,104)=1;
//r(eqn)=r(eqn)+Qv(gid(abs(jun_con(j,l+1)),jn_node))=250+Qv(52);
//??????????????///Therefore, equation becomes, 1*del_Qv(52)=250+Qv(52). But in my opinion, from equation 159 and 160, DBQ=Qv(52)+250 & d(DBQ)/d(Qv(52))=1. Using the relation, J*phi=-F(phi)(See 164 and 165), 1*del_Qv(52)=-(250+Qv(52))

/////For 2nd iteration: for j=2;
// eqn=eqn+1; 
//jun_inf(j,2)=jun_inf(1,2)=-99999;
// r(eqn)=0;
// l=1(i.e. for j=2 & l=1),
//if (abs(jun_con(j,l+1)))=abs(jun_con(2,2))=abs(4)=4 > eps_max;
// if (jun_con(j,l+1)=jun_con(2,2)=4>0)
// then jn_node=1;
//A(eqn,2*gid(abs(jun_con(j,l+1)),jn_node))=A(eqn,2*gid(4,1))=A(eqn,2*48)=A(eqn,96)=-1;
//r(eqn)=r(eqn)+Qv(gid(abs(jun_con(j,l+1)),jn_node))=250+Qv(gid(abs(jun_con(2,2)),1))=250+Qv(gid(4,1))=250-Qv(48);

//~~~~~~~~~~~Junction Energy condition~~~~~~~~

//    jun_con= [1 -4 0 0
//              3 4 -3 -2
//              3 -1 2 3];
for j=1:junn
    
    if(jun_inf(j,1)<>-99999) then
        eqn=eqn+1;
        if(jun_con(j,2)>0)then jn_nodel=1; end
        if(jun_con(j,2)<0)then jn_nodel=mnode(abs(jun_con(j,2))); end
        A(eqn, 2*gid(abs(jun_con(j,2)),jn_nodel)-1)=1;
        r(eqn)=yv(gid(abs(jun_con(j,2)),jn_nodel))-jun_inf(j,1);
        r(eqn)=-r(eqn); 
    end
    
    if(jun_con(j,2)>1) then
        for l=1:jun_con(j,1)-1
        eqn=eqn+1;
        if(jun_con(j,2)>0)then jn_nodel=1; end
        if(jun_con(j,2)<0)then jn_nodel=mnode(abs(jun_con(j,2))); end
        A(eqn, 2*gid(abs(jun_con(j,2)),jn_nodel)-1)=1;
        r(eqn)=yv(gid(abs(jun_con(j,2)),jn_nodel));
        
        if(jun_con(j,l+2)>0)then jn_node2=1; end
        if(jun_con(j,l+2)<0)then jn_node2=mnode(abs(jun_con(j,l+2))); end
        A(eqn, 2*gid(abs(jun_con(j,l+2)),jn_nodel)-1)=-1;
        r(eqn)=r(eqn)-yv(gid(abs(jun_con(j,l+2)),jn_node2));
        r(eqn)=-r(eqn); 
    end
end
end
//~~~~~~~~~Delta gv~~~~~~~
delyQ=A\r;
for i=1:2*sum(mnode)
    gv(i)=gv(i)+delyQ(i);
    rmse=rmse+delyQ(i)^2;
end
//~~~~~Update Value~~~~~
for l=1:chln
    for i=1:mnode(l)
        yv(gid(l,i))=gv(2*gid(l,i)-1);
        Qv(gid(l,i))=gv(2*gid(l,i));
    end
end
rmse=sqrt(rmse/sum(mnode));
count=count+1;
disp([count rmse])
end

//Print output
for l=1:chln
    disp(['channel number:' string(l)])
    disp('Section  Distance(m)  Depth(m)  Discharge(m^3/s)')
    for i=1:mnode(l)
    disp([i (i-1)*delta_x(l) yv(gid(l,i)) Qv(gid(l,i))])
end
end

////Wrong output!!!!!!!!!!!










































































