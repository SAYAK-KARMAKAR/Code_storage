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

function Cliv=Cli(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
//"y1o,Q1o,y2o,Q2o" repesents previous time level values.
    term1=areav(y2,B,m1,m2)-areav(y2o,B,m1,m2);
    term2=areav(y1,B,m1,m2)-areav(y1o,B,m1,m2);
    term3=Q2-Q1;
    term4=Q2o-Q1o;
    Cliv=(psi/delta_t)*term1+((1-psi)/delta_t)*term2+(theta/delta_x)*term3+((1-theta)/delta_x)*term4;
//See equation (194). As there is no extraction or injection at any junction; so "q" term in equation (194) was taken as zero.
endfunction

//~~~~~~~~~~Continuity Functions~~~~~~~~~~~~~~~~~~~
function dCdyiv=dCdyi(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
    dCdyiv=((1-psi)/delta_t)*dareav(y1,B,m1,m2);
    //See equation 195(.1)
endfunction

function dCdyip1v=dCdyip1(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
    dCdyip1v=((psi)/delta_t)*dareav(y2,B,m1,m2);
//See equation 195(.3)
endfunction

function dCdQiv=dCdQi(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
dCdQiv=-theta/delta_x;
//See equation 195(.2)    
endfunction

function dCdQip1v=dCdQip1(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
 dCdQip1v=theta/delta_x; 
//See equation 195(.4)    
endfunction

//~~~~~~~~~~~~~~Momentum Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function Mliv=Mli(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
    Av1=areav(y1,B,m1,m2);
    Av2=areav(y2,B,m1,m2);
    Av1o=areav(y1o,B,m1,m2);
    Av2o=areav(y2o,B,m1,m2);
//Av1,Av2 are corresponding to higher time level and  Av1o,Av2o are corresponding to Previous time level.   
    
    Rh1=HRv(y1,B,m1,m2);    
    Rh2=HRv(y2,B,m1,m2);
    Rh1o=HRv(y1o,B,m1,m2);
    Rh2o=HRv(y2o,B,m1,m2);
    
    term11=(Q2/Av2-Q2o/Av2o);
    term11=(Q1/Av1-Q1o/Av1o);
    term21=((alpha2/2)*Q2^2*Av2^(-2))-((alpha1/2)*Q1^2*Av1^(-2));
    term22=((alpha2/2)*Q2o^2*Av2o^(-2))-((alpha1/2)*Q1o^2*Av1o^(-2));
    term31=(y2+zv2)-(y1+zv1);
    term32=(y2o+zv2)-(y1o+zv1);
    term41=nm^2*Q2*abs(Q2)*Rh2^(-4/3)*Av2^(-2);
    term42=nm^2*Q1*abs(Q1)*Rh1^(-4/3)*Av1^(-2);
    term43=nm^2*Q2o*abs(Q2o)*Rh2o^(-4/3)*Av2o^(-2);
    term44=nm^2*Q1o*abs(Q1o)*Rh1o^(-4/3)*Av1o^(-2);
    
    Mliv=(psi/delta_t)*term11+((1-psi)/delta_t)*term12+(theta/delta_x)*term21+((1-theta)/delta_x)*term22+(theta*g/delta_x)*term31+((1-theta)*g/delta_x)*term32+theta*psi*g*term41+theta*(1-psi)*g*term42+(1-theta)*psi*g*term43+(1-theta)*(1-psi)*g*term44;
//See equation 196.
endfunction

function dMdyiv=dMdyi(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
    Av1=areav(y1,B,m1,m2);
    Rh1=HRv(y1,B,m1,m2); 
    dAv1=dareav(y1,B,m1,m2);
    dRh1=dHRv(y1,B,m1,m2);
    
    term1=(Q1/Av1^2)*dAv1;
    term2=(Q1^2/Av1^3)*dAv1;
    term3=theta*g/delta_x;
    term41=2*Q1*abs(Q1)*dAv1*Rh1^(-4/3)*Av1^(-3);
    term42=(4/3)*Q1*abs(Q1)*dRh1*Rh1^(-7/3)*Av1^(-2);
    dMdyiv=-((1-psi)/delta_t)*term1+(theta*alpha1/delta_x)*term2-term3-theta*(1-psi)*g*nm^2*(term41+term42);
//See equation 197(.1)
endfunction

function dMdyip1v=dMdyip1(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
    Av2=areav(y2,B,m1,m2);
    Rh2=HRv(y2,B,m1,m2); 
    dAv2=dareav(y2,B,m1,m2);
    dRh2=dHRv(y2,B,m1,m2);
    
    term1=(Q2/Av2^2)*dAv2;
    term2=(Q2^2/Av2^3)*dAv2;
    term3=theta*g/delta_x;
    term41=2*Q2*abs(Q2)*dAv2*Rh2^(-4/3)*Av2^(-3);
    term42=(4/3)*Q2*abs(Q2)*dRh2*Rh2^(-7/3)*Av2^(-2);
    dMdyip1v=-(psi/delta_t)*term1-(theta*alpha2/delta_x)*term2+term3-theta*psi*g*nm^2*(term41+term42);
//See equation 197(.3)
endfunction

function dMdQiv=dMdQi(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
    Av1=areav(y1,B,m1,m2);
    Rh1=HRv(y1,B,m1,m2); 
    term1=Av1^(-1);
    term2=Q1/Av1^2;
    term3=abs(Q1)*Rh1^(-4/3)*Av1^(-2);
    dMdQiv=((1-psi)/delta_t)*term1-(theta*alpha1/delta_x)*term2+2*theta*(1-psi)*g*nm^2*term3;
//See equation 197(.2)
endfunction

function dMdQip1v=dMdQip1(y1,Q1,y2,Q2,y1o,Q1o,y2o,Q2o,zv1,zv2,theta,psi,delta_t,delta_x,alpha1,alpha2,B,m1,m2,nm)
    Av2=areav(y2,B,m1,m2);
    Rh2=HRv(y2,B,m1,m2); 
    term1=Av2^(-1);
    term2=Q2/Av2^2;
    term3=abs(Q2)*Rh2^(-4/3)*Av2^(-2);
    dMdQip1v=(psi/delta_t)*term1+(theta*alpha2/delta_x)*term2+2*theta*psi*g*nm^2*term3;
//See equation 197(.4)
endfunction

//~~~~~~~~Boundary values Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function bndv=bndcon(typ,jnum,tv)
//'typ' means type. if typ=1, means we specify depth. If typ=2, we specify the discharge value. "jnum" specify the junction number(See 5 no-page-49)."tv" means time value.

if (typ==1 & jnum==3) then
    bndv=1.43;
// Depth at downstream has a fixed value.
end

if (typ==2 & jnum==1) then
    if (tv<2000) then
        bndv=50+100/2000*tv;
    end
    if (tv>=2000) then
        bndv=150-(100/2000)*(tv-2000);
    end
    if (tv>4000) then
        bndv=50;
    end     
end

if (typ==2 & jnum==2) then
    if (tv<2000) then
        bndv=50+100/2000*tv;
    end
    if (tv>=2000) then
        bndv=150-(100/2000)*(tv-2000);
    end
    if (tv>4000) then
        bndv=50;
    end 
end
// At junction number 1 and 2, there are varying discharge boundary conditions with time.
endfunction
 

// Channel reach: Start + End -
//~~~~~~~~~~~~~~Given Data~~~~~~~~~
g=9.81;//m.s^(-2)
global('g')
eps_max=1e-6;
t_max=25001;//s
delta_t=250;//s
theta=0.5;
// To get Intermediate value between 'n' and 'n+1'th time level.
psi=0.5;
// To get Intermediate value between 'i' and 'i+1'th node.
junnn=4;
bjn=3;
chln=3;
//~~~~Chl | Length | Width | m1 | m2 | Segment | n | S0 | JN1 | JN2 |
chl_inf=[1 5000  50 0 0 500 0.025 0.0002 1 4
         2 5000  50 0 0 500 0.025 0.0002 2 4
         3 5000 100 0 0 500 0.025 0.0002 4 3];
//~~~~~~Specified flow depth | Specified Discharge | Bed elevation~~~~~~
jun_inf=[-99999 2 2
         -99999 2 2
         1 -99999 0
         -99999 -99999 1];
//"1" is symbolized as specified depth."2" is symbolized as specified discharge."-99999" means no value specified.
jun_con=[1 1 0 0
         1 2 0 0
         1 -3 0 0
         3 3 -1 -2];
//In 'jun_con' matrix, 1st column represents number of channels connected to that junction.
//Positive sign means 1st section of that channel is connected to that junction. Negative sign means (Nl+1)th section of that channel is connected to that junction. 
alpha=ones(chln,1);

//Derived informations
Lx=chl_inf(1:chln,2);
//We know, 'chln' means channel number=4.
B=chl_inf(1:chln,3);
m1=chl_inf(1:chln,4);
m2=chl_inf(1:chln,5);
delta_x=chl_inf(1:chln,6);
nm=chl_inf(1:chln,7);
S0=chl_inf(1:chln,8);

mnode=Lx./delta_x+1;
//mnode is total number of sections for a particular channel. Here, mnode is calculated at point by point basis.

//~~~~~~~~z values~~~~~~~~~~
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

//~~~~~~~Problem Dependent Parameters~~~~~~~
yv=zeros(sum(mnode),1);
Qv=zeros(sum(mnode),1);
yo=zeros(sum(mnode),1);
Qo=zeros(sum(mnode),1);
// yv, Qv are values corresponding to future time level and yo, Qo are values corresponding to old time level.

//~~~General identification matrix~~~
idv=0;
for l=1:chln
    for i=1:mnode(l)
        idv=idv+1;
        gid(l,i)=idv;
//Number of rows in 'gid' matrix= no. of channels=4 and number of columns= maximum no. of sections in a channel reach=max(11,11,11)=11. So, 'gid' is a 3*11 matrix.
    end
end

//~~~~~~~Initial values~~~~~~~~~
for l=1:chln
    for i=mnode(l)
        if (l==1 & l==2) then
            yo(gid(l,i))=1.4300;
            Qo(gid(l,i))=50.00;
// Initial depth and discharge values for channel 1 and 2.
        else
            yo(gid(l,i))=1.4300;
            Qo(gid(l,i))=100.00;
// Initial depth and discharge values for channel 3.
        end
    end
end

yv=yo;
Qv=Qo;
//At first, yv and Qv vector stores initial values.
gv=zeros(2*sum(mnode),1);

for l=1:chln
     for i=1:mnode(l)
     gv(2*gid(l,i)-1)=yv(gid(l,i));
     gv(2*gid(l,i))=Qv(gid(l,i));
//Initial Values of 'yv' and 'Qv' matries are stored in 'gv' or general variable matrix togetherly.
     end
end
//~~~~~MY input~~~~~~~
//timev=zeros((t_max-1)/delta_t,1);
//ysim=zeros((t_max-1)/delta_t,1);
//Qsim=zeros((t_max-1)/delta_t,1);
//~~~~~~~~~~~~~~~~~~ Time loop~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tv=0;//zero time level.
tcount=0;
while tv<t_max
    tv=tv+delta_t;
    A=zeros(2*sum(mnode),2*sum(mnode));
//Specifying the size of the jacobian matrix.
    r=zeros(2*sum(mnode),1);
    disp(['Time in seconds:' string(tv)])
    count=0;
    rmse=1;
    
    //~~~~~~~~Space Loop~~~~~~~~~~
    while rmse>eps_max
        rmse=0;
//????? In line 286, rmse=1 and here rmse=0. What does it mean???
        eqn=0;
//Equation number.

 //~~~~Equations corresponding to segments (2N1+2N2+2N3)~~~~
    for l=1:chln
        for i=1:mnode(l)-1
            //~~~~~~~Continuity~~~~~~~
            eqn=eqn+1;
            A(eqn,2*gid(l,i)-1)=dCdyi(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(C_{l,i})/del(y_{l,i}).See page 52 & equation 195(.1). Everything written within the brackets are variables corresponding to i-th segment of l-th channel (i.e.(l,i)).
         A(eqn,2*gid(l,i))=dCdQi(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(C_{l,i})/del(Q_{l,i}).See equation 195(.2).
          A(eqn,2*gid(l,i+1)-1)=dCdyip1(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(C_{l,i})/del(y_{l,i+1}).See equation 195(.3).
          A(eqn,2*gid(l,i+1))=dCdQip1(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(C_{l,i})/del(Q_{l,i+1}).See equation 195(.4).
          
          r(eqn)=-Cli(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
            //~~~~~~~Momentum~~~~~~~
            eqn=eqn+1;
            A(eqn,2*gid(l,i)-1)=dMdyi(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(M_{l,i})/del(y_{l,i}).See page 52 & equation 197(.1). Everything written within the brackets are variables corresponding to i-th segment of l-th channel (i.e.(l,i)).
         A(eqn,2*gid(l,i))=dMdQi(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(M_{l,i})/del(Q_{l,i}).See equation 197(.2).
          A(eqn,2*gid(l,i+1)-1)=dMdyip1(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(M_{l,i})/del(y_{l,i+1}).See equation 197(.3).
          A(eqn,2*gid(l,i+1))=dMdQip1(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
//del(M_{l,i})/del(Q_{l,i+1}).See equation 197(.4).
          
          r(eqn)=-Mli(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),yo(gid(l,i)),Qo(gid(l,i)),yo(gid(l,i+1)),Qo(gid(l,i+1)),zv(l,i),zv(l,i+1),theta,psi,delta_t,delta_x(l),alpha(l),B(l),m1(l),m2(l),nm(l));
          end
end

//~~~~~~Junction Continuity condition~~~~~~~
for j=1:junn
    dcond=0; //Initially zero,=1 if discharge condition is there.
    if (jun_inf(j,2)==2) then
// At the 2nd column of jun_inf matrix, discharge values of the junctions are shown. "2" number symbolizes discharge and -99999 represents no specified value.
        eqn=eqn+1;
        r(eqn)=bndcon(2,j,tv);
//We know "bndv=bndcon(typ,jnum,tv)". For example bndcon(2,1,0) means discharge boundary condition of 1st junction at time=0s.
        dcond=1;
    else
        if(j>bjn) then
//"bjn" means number of boundary junctions(i.e=3). Those are numbered as 1,2 and 3. So, this condition implies the internal junctions. 
            eqn=eqn+1;
            r(eqn)=0;
            dcond=1;
        end
    end
    
    if (dcond==1) then
        for l=1:jun_con(j,1)
            if(abs(jun_con(j,l+1))>eps_max) then
                if (jun_con(j,l+1)>0) then 
                    jun_node=1;
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
end

//~~~~~~~~~Junction Energy Condition~~~~~~~~~~
for j=1:junn
    
    if(jun_inf(j,1)==1) then
//1st column of "jun_inf" matrix represents the specified depth at boundary junctions (if specified,represented as 1. If unspecified, represented as -99999).
        eqn=eqn+1;
        if(jun_con(j,2)>0)then jn_nodel=1; end
        if(jun_con(j,2)<0)then jn_nodel=mnode(abs(jun_con(j,2))); end
//Above two 'if' statements check whether the junction is connected to 1st section or the last section of the channel. 
        A(eqn, 2*gid(abs(jun_con(j,2)),jn_nodel)-1)=1;
        r(eqn)=yv(gid(abs(jun_con(j,2)),jn_nodel))-bndcon(1,j,tv);
        r(eqn)=-r(eqn); 
    end
    
//Explanation of the above loop:
//Let us start with j=3. Because Only for that We have specified depth condition.
//for j=3;
//if(jun_inf(j,1)== jun_inf(3,1)==1) then 
//       eqn=eqn+1; 
//       if(jun_con(j,2)=jun_con(3,2)=-3<0)
//            then jn_nodel=mnode(abs(jun_con(j,2)))=mnode(abs(jun_con(3,2)))=mnode(abs(-3))=mnode(3)=11; 
//        end    
//        A(eqn, 2*gid(abs(jun_con(j,2)),jn_nodel)-1)=A(eqn, 2*gid(abs(jun_con(3,2)),jn_nodel)-1)=A(eqn, 2*gid(3,11)-1)=A(eqn,2*33-1)=A(eqn,65)=1;
//         r(eqn)=yv(gid(abs(jun_con(j,2)),jn_nodel))-bndcon(1,j,tv)=yv(gid(abs(jun_con(3,2)),jn_nodel))-bndcon(1,j,tv)=yv(gid(3,11))-bnd(1,3,0)=yv(33)-1.43;
//          r(eqn)=-r(eqn)=1.43-yv(33)
//So, as per this loop, d/s equation becomes, 1.delta_y(33)=1.43-yv(33).

//?????????? But in my opinion,Mli=yd-y(33), dMdyi=-1. So, acoording to the rule dMdyi*delta_y(33)=-Mli => (-1).delta_y(33)=-(yd-y(33)) => "delta_y(33)=y(33)-yd".    this is different!

    if(jun_con(j,1)>1) then
//It means, if more than one channel is connected to that junction.
        for l=1:jun_con(j,1)-1
//It runs for (number of channels connected to that junction-1) times. Because, one channel is already calculated for the channel represented by "jun_con(j,2)" elements in previous loop.
        eqn=eqn+1;
//"
        if(jun_con(j,2)>0)then jn_nodel=1; end
        if(jun_con(j,2)<0)then jn_nodel=mnode(abs(jun_con(j,2))); end
        A(eqn, 2*gid(abs(jun_con(j,2)),jn_nodel) -1)=1;
        r(eqn)=yv(gid(abs(jun_con(j,2)),jn_nodel));
//" 
//?????????Why "  "  part is required? I think, it has been already calculated. Ans: No, junction with multiple channel is seperately calculated. This part is exclusively for junction 4.    
        if(jun_con(j,l+2)>0)then jn_node2=1; end
// Minimum and maximum values of l are:(1 & jun_con(j,1)-1). So, it will run for 1+2=3rd column to last non-zero column (depends on how many channels are connected to that particular junction) of "jun_con" matrix.
        if(jun_con(j,l+2)<0)then jn_node2=mnode(abs(jun_con(j,l+2))); end
        A(eqn, 2*gid(abs(jun_con(j,l+2)),jn_node2)-1)=-1;
        r(eqn)=r(eqn)-yv(gid(abs(jun_con(j,l+2)),jn_node2));
        r(eqn)=-r(eqn); 
    end
end
end

//~~~~~del y Q~~~~~~~~
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
//We should update those values so that we can use those values for the next iteration.
    end
end
rmse=sqrt(rmse/sum(mnode));
count=count+1;
//disp([count rmse])
end

tcount=tcount+1;
//~~~~~~Update values for time n<--n+1~~~~~~
for l=1:chln
for i=1:mnode(l)
    yo(gid(l,i))=yv(gid(l,i));
    Qo(gid(l,i))=Qv(gid(l,i));
end
end

timev(tcount)=tv;
ysim(tcount)=yv(31);
Qsim(tcount)=Qv(31);

end

//Print output
plot(timev,ysim)
xtitle('Flow depth at 4000 m in channel 3','Time(s)', 'Flow Depth(m)');
figure
plot(timev,Qsim)
xtitle('Discharge at 4000 m in channel 3','Time(s)', 'Discharge(m^3/s)');
























































































































