clc
clear
//Lecture 39. variables: y and Q Both.
//~~~~~~Given Data~~~~~~~
junc=1;
chl=2;
QI=20;//m^3/s
S0=[0.0004 0.0008];
n=[0.010 0.015];
B=15;//m
g=9.81;//m/s^2
Lx=[100 100];
yd=0.6;//m
//Depth at the downstream section.
mnode=[101 101];
eps_max=1e-6;
global('B','g')
//here two channels are there in series. So, C1,C2,S0,n values are different for the channels. Also Q is considered as variable, no constant dischargefor this problem.
juni=[1 2 101 1]; 
//I.T-->Coordinate of the junction =(channel number, node number)=(1,101)=(2,1).

//~~~~~~~~Problem Dependent Parameters~~~~~~~~~~
alpha=[1,1];
yv=yd*ones(sum(mnode),1);
Qv=QI*ones(sum(mnode),1);
//We will initialize the problem with yd and QI values of depth and discharge.Both are (202*1) matrices.
gv=zeros(2*sum(mnode),1);
//"gv" is the general variable with both y and Q.

//~~~~~~General Identification matrix~~~~~
 idv=0;
 for l=1:chl
     for i=1:mnode(l)
         idv=idv+1
         gid(l,i)=idv
//Previously, we used local numbering of a section i.e.(l,i).Now, this local numberings are converted to the global numbering and stored in 'gid' matrix. Eg. for l=chl=2, for i=mnode(2)=101, gid(2,101)=202.
     end
 end
 for l=1:chl
     for i=1:mnode(l)
     gv(2*gid(l,i)-1)=yv(gid(l,i))
     gv(2*gid(l,i))=Qv(gid(l,i))
     //Here, yv and Qv vectors store depth and dischage informations of the sections (size 202*1 for each). Where, 'gv'(size 404*1) is the "general variable" vector which collects data from 'yv' and 'Qv' vectors. 
     //Here, equation looks like, "J*del_phi=r".Whare, J=Jacobian matrix, 
    //del_phi={del_y1; del_Q1; del_y2; del_Q2;......},
    //r=[-UB_{1,1}; -C_{1,1}; -M_{1,1}; -C_{1,2}; -M_{1,2}; .........;-C_{chl,mnode(chl)}; -M_{chl,mnode(chl)};-DB_{chl,mnode(chl)}]
     end
 end

for l=1:chl
    delta_x(l)=Lx(l)/(mnode(l)-1);
    D1(l)=alpha(l)/(2*g);
    D2(l)=(1/2)*n(l)^2*delta_x(l);
//These three are channel reach dependent parameters.
end
mc=sum(mnode);
for l=chl:-1:1
    for i=mnode(l):-1:1
        if (l==chl & i==mnode(l)) then
            zv(mc)=0;
//This 1st 'if' loop evaluates the bed elevation for end section of the last channel reach.
        end
        if(l<>chl & i==mnode(l)) then
            mc=mc-1;
            zv(mc)=zv(mc+1);
//This 'if' loop evaluates the bed elevation for end section of other channel reaches.
        end
        if(i<>mnode(l)) then
            mc=mc-1;
            zv(mc)=zv(mc+1)+S0(l)*delta_x(l);
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
function Mliv=Mli(y1,Q1,y2,Q2,S0,delta_x,D1,D2)

    Mliv=(y2-y1)-S0*delta_x+D1*(Q2*abs(Q2)*areav(y2)^(-2)-Q1*abs(Q1)*areav(y1)^(-2))+D2*(Q2*abs(Q2)*HRv(y2)^(-4/3)*areav(y2)^(-2)+Q1*abs(Q1)*HRv(y1)^(-4/3)*areav(y1)^(-2));
 //See Equation (146).
 //Here, coefficients are D1 and D2, not C1 and C2. C1 and C2 was considered for the code "steady_1D_channel_series", Where 'discharge Q' was not considered as a variable. 
 //////////////Why Q|Q| in Kineric head term?????  
endfunction

function dMdyiv=dMdyi(y,Q,D1,D2)
    term1=(2*Q^2/areav(y)^(3))*dareav(y);
    term2=2*Q^2*areav(y)^(-3)*HRv(y)^(-4/3)*dareav(y);
    term3=(4/3)*Q^2*areav(y)^(-2)*HRv(y)^(-7/3)*dHRv(y);
    dMdyiv=-1+D1*term1-D2*(term2+term3);
    //From equation (147.1).Evaluates dM_{l,i}/dy_{l,i}.
endfunction

function dMdyip1v=dMdyi(y,Q,D1,D2)
    term1=(2*Q^2/area(y)^(3))*dareav(y);
    term2=2*Q^2*areav(y)^(-3)*HRv(y)^(-4/3)*dareav(y);
    term3=(4/3)*Q^2*areav(y)^(-2)*HRv(y)^(-7/3)*dHRv(y);
    dMdyiv=-1-D1*term1-D2*(term2+term3);
    //From equation (147.3). Evaluates dM_{l,i}/dy_{l,i+1}.
endfunction
    
function dMdQip1v=dMdQip1(y,Q,D1,D2)
    term1=2*Q*areav(y)^(-3);
    term2=2*Q*areav(y)^(-2)*HRv(y)^(-4/3);
    dMdQip1v=D1*term1+D2*term2;
    //From equation (147.4).
endfunction

function dMdQiv=dMdQi(y,Q,D1,D2)
    term1=2*Q*areav(y)^(-3);
    term2=2*Q*areav(y)^(-2)*HRv(y)^(-4/3);
    dMdQip1v=-D1*term1+D2*term2;
    //From equation (147.2).
endfunction

A=zeros(2*sum(mnode),2*sum(mnode));
r=zeros(2*sum(mnode),1);
// We have total 101*2=202 nodes in two channels. So, we have 202*2=404 variables in total (considering y & Q).
count=0;
rmse=1;

//~~~~~~~~~~~~~~~~~~~~~~~~~~Space Loop~~~~~~~~~~~~~~~~~~~~~~~~~
while rmse>eps_max
    rmse=0;
    eqn=1; 
//Equation number.
//~~~Upstream Boundary~~~~
    A(1,2)=1;
    r(1)=-(Qv(1)-QI);
    //See (5no-P 133). Here, UB_{1,1}=Q_{1,1}-Q_u. A(1,1)=del(UB)/del_(y_{1,1})=0. So, this satisfy "A*phi=r", means "Jacobian*del_(y and Q)=-(C_{l,i} and M_{l,i})". See (5no-P 33).
//~~~~~~~~~~~Equations Corresponding to segments~~~~~~~~~~
for l=1:chl
    for i=1:mnode(l)-1
// Because, we get one continuity and one momentum equation corresponding to each segment. and for lth channel reach, number of segments=mnode(l)-1. Here, (l,i) means the segment number which have (l,i)th and (l,i+1)th sections.
        eqn=eqn+1;
//~~~~~~Jacobians for Continuity eqn~~~~~~~
        A(eqn,2*gid(l,i)-1)=0;
//del(C_{l,i})/del(y_{l,i}).See equation (156).
        A(eqn,2*gid(l,i))=-1;
//del(C_{l,i})/del(Q_{l,i}).
        A(eqn,2*gid(l,i+1)-1)=0;
//del(C_{l,i})/del(y_{l,i+1}).
        A(eqn,2*gid(l,i+1))=1;
//del(C_{l,i})/del(Q_{l,i+1}).
        r(eqn)=0;
// Value of each continuity functions (i.e, C_{l,i}) is zero! See equation (143) & (155.1).
        eqn=eqn+1;
//This 'increment of counter' is for 'momentum equation' of the same segment.
        
//~~~~~~~Jacobians for Momentum eqn~~~~~~~~
        A(eqn,2*gid(l,i)-1)=dMdyi(yv(gid(l,i)),Qv(gid(l,i)),D1(l),D2(l));
//del(M_{l,i})/del(y_{l,i}).See page 37 & equation (147.1).
         A(eqn,2*gid(l,i))=dMdQi(yv(gid(l,i)),Qv(gid(l,i)),D1(l),D2(l));
//del(M_{l,i})/del(Q_{l,i}).See equation (147.2).
          A(eqn,2*gid(l,i+1)-1)=dMdyip1(yv(gid(l,i+1)),Qv(gid(l,i+1)),D1(l),D2(l));
//del(M_{l,i})/del(y_{l,i+1}).See equation (147.3).
          A(eqn,2*gid(l,i+1))=dMdQip1(yv(gid(l,i+1)),Qv(gid(l,i+1)),D1(l),D2(l));
//del(M_{l,i})/del(Q_{l,i+1}).See equation (147.4).
          
          r(eqn)=-Mli(yv(gid(l,i)),Qv(gid(l,i)),yv(gid(l,i+1)),Qv(gid(l,i+1)),S0(l),delta_x(l),D1(l),D2(l));
          end
end
// Explanation of 'equations corresponding to segments' of 'space loop': 
//1st iteration (When l=1 & i=1):       Continuity part:eqn=1+1=2 (i.e. 2nd row of Jacobian matrix). gid(1,1)=1 and gid(1,2)=2. Therefore, A(2,1)=0; A(2,2)=-1; A(2,3)=0; A(2,4)=1.        Momentum part: eqn=2+1=3 (i.e. 3rd row of Jacobian matrix). gid(1,1)=1. To evaluate A(3,1), A(3,2), A(3,3)and A(3,4) equation (147)'s are used.
//2nd iteration(When l=1 & i=2): then eqn=4. gid(1,2)=2 and gid(1,3)=3. Therefore, A(4,3)=0; A(4,4)=-1; A(4,5)=0; A(4,6)=1....and so on...

//~~~~~~~~~Junction Condition~~~~~~
///////Doubt: If junction condition executes after all segments, is the position of junction equation not already occupied???

//~~~~Junction Continuity~~~
eqn=eqn+1;
A(eqn,2*gid(juni(1),juni(3)))=-1;
A(eqn,2*gid(juni(2),juni(4)))=1;
r(eqn)=-(Qv(gid(juni(2),juni(4)))-Qv(gid(juni(1),juni(3))));
////I.T=>[eqn=203 should be here(Because, previously 202 eqns are there i.e. 1 u/s condition and 200 equations for 1st 100 segments and 1 junction energy condition (because, delta_y comes before delta_Q.).
//A(202,2*gid(1,101))=A(202,2*101)=A(202,202)=-1 ;
//A(202,gid(2,1))=A(202,2*102)=A(202,204)=1 and
//r(202)=-(Qv(102)-Qv(101))=-(Qv_{2,1}-Qv_{1,101})
//So, equation becomes, -delta_Qv(101)+delta_Qv(102)=-(Qv(102)-Qv(101)). ****Remember, delta_Qv(101), delta_Qv(102) these are difference of current and previous iteration values from Newton-Raphson's method, not the difference between two node values.]

//~~~~Junction Energy~~~
eqn=eqn+1;
A(eqn,2*gid(juni(1),juni(3))-1)=1;
A(eqn,2*gid(juni(2),juni(4))-1)=-1;
r(eqn)=-(yv(gid(juni(1),juni(3)))-yv(gid(juni(2),juni(4))));
/////I.T--> eqn=202 should be here. Because, junction energy comes before junction continuity (i.e, delta_y comes before delta_Q). But, in code it gets more counter(as, eqn=eqn+1 used). Is it right?????

////I.T--> [We have total 101+101=202 sections=>202*2=404 variables. We have 1 upstream discharge condition; 100+100=200 segments=>200*2=400 equations(Continuity & momentum); 2 junction conditions (Continuity & momentum); 1 downstrem depth condition i.e. 1+400+2+1=404 equations].

//~~~~~Downstream Boundary~~~~~~~~~
eqn=eqn+1;
A(eqn,2*gid(2,mnode(2))-1)=-1;
// d/s depth boundary, coff corresponding to delta_y(202).
r(eqn)=-(yd-yv(gid(2,mnode(2))));
//So, d/s becomes, -delta_y(202)=-(yd-yv(202))

delyQ=A\r;
for i=1:2*sum(mnode) 
//Because, 2*sum(mnode)=2*202=404 is total number of variables.
    gv(i)=gv(i)+delyQ(i); 
//All the variables are being updated in each iteration.
    rmse=rmse+delQ(i)^2;
end

//Initial Value
for l=1:chl
    for i=1:mnode(l)
        yv(gid(l,i))=gv(2*gid(l,i)-1);
        Qv(gid(l,i))=gv(2*gid(l,i));
    end
end
    rmse=sqrt(rmse/sum(mnode));
    count=count+1;
    disp([count rmse])
end
//figure
plot(xv,yv+zv,"-r")
plot([Lx(1) Lx(1)+Lx(2)],[zv(mnode(1)+1) zv(sum(mnode))],'m-')
plot([0,Lx(1)],[zv(1),zv(mnode(1))],'m-')
xtitle('steady channel flow', 'x axis','Flow depth')




































































