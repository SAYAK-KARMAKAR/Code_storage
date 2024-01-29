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
yd=5; //m
Qd=250; //m^3/s
Qu=250; //m^3/s
eps_max=1e-6;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
junn=6;//number of junctions
bjn=2;//Out of 6 junctions, 2 are boundary junctions.
chln=8;//number of channels
//~~~~Chl | Length | Width | m1 | m2 | Segment | n | S0 | JN1 | JN2 |~~~~~~~~
chl_inf=[1  200 30 0 0 50 0.0130 0.0005 1 3
         2  200 40 0 0 50 0.0130 0.0005 1 6
         3  200 20 0 0 50 0.0120 0.0005 3 4
         4  100 20 0 0 25 0.0140 0.0005 3 5
         5  100 20 0 0 25 0.0130 0.0005 5 4
         6  100 25 0 0 25 0.0130 0.0005 6 5
         7  100 30 0 0 25 0.0140 0.0005 4 2
         8  300 50 0 0 75 0.0140 0.0005 6 2];
//~~~~~~Specified flow depth | Specified Discharge | Bed elevation~~~~~~
jun_inf=[-99999 Qu 0.25
         yd -Qd 0
         -99999 -99999 0.15
         -99999 -99999 0.05
         -99999 -99999 0.1
         -99999 -99999 0.15];
//Upstream junction is numbered as 1 & downstream junction as 2.
jun_con= [2 1 2 0
          2 -7 -8 0
          3 -1 3 4
          3 -3 -5 7
          3 -4 -6 5
          3 -2 6 8];
alpha=ones(chln,1);
//~~~~~~~Derived informations from 'chl_inf' matrix~~~~~~~~~~
Lx=chl_inf(1:chln,2);
//We know, 'chln' means channel number=8.
B=chl_inf(1:chln,3);
m1=chl_inf(1:chln,4);
m2=chl_inf(1:chln,5);
delta_x=chl_inf(1:chln,6);
n=chl_inf(1:chln,7);
S0=chl_inf(1:chln,8);

mnode=Lx./delta_x+1;
//Here, mnode is calculated at point by point basis.

//~~~~~~~~~~~z values~~~~~~~~~~~~~
for l=1:chln
    if (jun_inf(chl_inf(l,9),3)>jun_inf(chl_inf(l,10),3)) then
        fact=-1;
// 9th column of 'chl_inf' matrix represents the junction connected to the 1st section of that particular channel. 10th column represents the 'junction' connected to the last section of that particular channel. It checks, which junction is at higher elevation between this two...
    else
        fact=+1;
end
zv(l,1)=jun_inf(chl_inf(l,9),3);
for i=2:mnode(l)
    zv(l,i)=zv(l,i-1)+fact*S0(l)*delta_x(l);
//Determines the elevations of the intermediate sections of l-th channel.
    end
end
///Explanation of the 'z values' loop:
//1st iteration:
//for l=1:chln=1:8=> Take,l=1
//    if (jun_inf(chl_inf(1,9),3)=(jun_inf(1,3))=0.25 > jun_inf(chl_inf(1,10),3))=(jun_inf(3,3))=0.15 then
//        fact=-1;
//end
//zv(1,1)=jun_inf(chl_inf(1,9),3)=jun_inf(1,3)=0.25;
//for i=2:mnode(1)=2:5 (take i=2)
//    zv(1,2)=zv(1,i-1)+fact*S0(1)*delta_x(1)=zv(1,1)+(-1)*S0(1)*delta_x(1)=0.25-0.0005*50=0.225;
//    end
//end
//......and so on.....


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
//Number of rows in 'gid' matrix= no. of channels=8 and number of columns= maximum no. of sections in a channel reach=max(5,5,5,5,5,5,5,5)=5. So, 'gid' is a 8*5 matrix. Unlike the previous code (Steady_1D_channel_network_Without_reverse), we have no 'zero" entry in 'gid' matrix. So, we have total 40 sections.
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
    
    //~~~~~Equations corresponding to segments(2N1+2N2+2N3+2N4+2N5+2N6+2N7+2N8)~~~~~~
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
jeqn=0;

//~~~~~~Junction Continuity Condition~~~~~~~
for j=1:junn
    if (jun_inf(j,2)~=-99999) then
//Checks whether the junction has specified didschsrge value or not.then r(eqn)=That specified discharge value.
        eqn=eqn+1;
        jeqn=jeqn+1;
        r(eqn)=jun_inf(j,2);
    else
        if(j>bjn) then
//'bjn' means number of boundary junction and j>bjn means it is applicable for internal junctions where no specified depth or discharge condition exists.
            eqn=eqn+1;
            jeqn=jeqn+1;
            r(eqn)=0;
    end
end

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


//~~~~~~~~~~~Junction Energy condition~~~~~~~~
//    jun_con= [1 -4 0 0
//              3 4 -3 -2
//              3 -1 2 3];
for j=1:junn
    
    if(jun_inf(j,1)<>-99999) then
//1st column of "jun_inf" matrix represents the specified depth at boundary junctions (if unspecified, represented as -99999).
        eqn=eqn+1;
        if(jun_con(j,2)>0)then jn_nodel=1; end
        if(jun_con(j,2)<0)then jn_nodel=mnode(abs(jun_con(j,2))); end
        A(eqn, 2*gid(abs(jun_con(j,2)),jn_nodel)-1)=1;
        r(eqn)=yv(gid(abs(jun_con(j,2)),jn_nodel))-jun_inf(j,1);
        r(eqn)=-r(eqn); 
    end
    
    if(jun_con(j,1)>1) then
//It means, if more than one channel is connected to that junction.
        for l=1:jun_con(j,1)-1
//It runs for (number of channels connected to that junction-1) times. Because, one channel is already calculated for the channel represented by "jun_con(j,2)" elements in previous loop.
        eqn=eqn+1;
//"
        if(jun_con(j,2)>0)then jn_nodel=1; end
        if(jun_con(j,2)<0)then jn_nodel=mnode(abs(jun_con(j,2))); end
        A(eqn, 2*gid(abs(jun_con(j,2)),jn_nodel)-1)=1;
        r(eqn)=yv(gid(abs(jun_con(j,2)),jn_nodel));
//" 
//?????????Why "  "  part is required? I think, it has been already calculated.     
        if(jun_con(j,l+2)>0)then jn_node2=1; end
// Minimum and maximum values of l are:(1 & jun_con(j,1)-1). So, it will run for 1+2=3rd column to last non-zero column (depends on how many channels are connected to that particular junction) of "jun_con" matrix.
        if(jun_con(j,l+2)<0)then jn_node2=mnode(abs(jun_con(j,l+2))); end
        A(eqn, 2*gid(abs(jun_con(j,l+2)),jn_node2)-1)=-1;
        r(eqn)=r(eqn)-yv(gid(abs(jun_con(j,l+2)),jn_node2));
        r(eqn)=-r(eqn); 
    end
end
end

//~~~~Del y Q~~~
delyQ=inv(A'*A)*(A'*r);
// Here, we have 80 unknowns. But we get 81 equations. (See '5no-P-43' to know why). So we get,
//    [A]_(81*80) * [delyQ]_(80*1)=[r]_(81*1)
// But, this matrix operation is not possible.So, multiplying both  side of the equation with 'transpose of A',
//    [A']_(80*81) * [A]_(81*80) * [delyQ]_(80*1)=[A']_(80*81) * [r]_(81*1)
//    [A'A]_(80*80) * [delyQ]_(80*1)=[A'r]_(80*1)
//     [delyQ]_(80*1)= inv([A'A]_(80*80)) * [A'r]_(80*1)


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
disp(eqn)
//Print output
for l=1:chln
    disp(['channel number:' string(l)])
    disp('Section  Distance(m)  Depth(m)  Discharge(m^3/s)')
    for i=1:mnode(l)
    disp([i (i-1)*delta_x(l) yv(gid(l,i)) Qv(gid(l,i))])
end
end



//Wrong Output!!!!!!!!!!!!!!!






































































