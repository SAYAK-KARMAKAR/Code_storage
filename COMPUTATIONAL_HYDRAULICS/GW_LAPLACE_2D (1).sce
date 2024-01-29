//GW LAPLACE 2D
clc
clear
//...problem dependent parameters...
mnode=31;
nnode=21;
Lx=300;
Ly=100;
hA=90;
hB=89;
hC=85;
hD=87;
//calculated values of parameters
delta_x=Lx/(mnode-1); 
delta_y=Ly/(nnode-1);
x=0:delta_x:Lx;
y=0:delta_y:Ly;
alphax=1/(delta_x)^2;
alphay=1/(delta_y)^2;
mnnode=mnode*nnode;
A=zeros(mnnode,mnnode);
r=zeros(mnnode,1);
for j=1:nnode
for i=1:mnode
l=i+(j-1)*mnode; //single index notation
//corner points are considered at specified boundary condition
//node A
if (i==1 & j==nnode) then
A(l,l)=1;
r(l)=hA;
end
//node B
if (i==1 &j==1) then
A(l,l)=1;
r(l)=hB;
end
//Node c
if (i==mnode & j==1)then
A(l,l)=1;
r(l)=hC;
end
//Node D
if (i==mnode & j==nnode)then
A(l,l)=1;
r(l)=hD;
end
//Interior points
if (i>1 & i<mnode) then
if(j>1 & j<nnode)then
    A(l,l-mnode)=alphay;
    A(l,l-1)=alphax;
    A(l,l)=-2*(alphax+alphay)
    A(l,l+1)=alphax;
    A(l,l+mnode)=alphay;
    r(l)=0;
end
end
//specified LBC
if (i==1) then
    if (j>1 & j<nnode)then
        A(l,l)=1.0;
        r(l)=hB+(hA-hB)*(j-1)*(delta_y/Ly);
    end
end
//specified RBC
if (i==mnode) then
    if (j>1 & j<nnode)then
        A(l,l)=1.0;
        r(l)=hC+(hD-hC)*(j-1)*(delta_y/Ly);
    end
end
//No flow BBC
if (j==1)then
    if (i>1 & i<mnode)then
        //3 point boundary
        A(l,l)=-3/(2*delta_y);
        A(l,l+mnode)=4/(2*delta_y);
        A(l,l+2*mnode)=-1/(2*delta_y);
        r(l)=0;
    end
end
//No flow TBC
if (j==nnode) then
    if (i>1 & i<mnode)then
        //3 point boundary
        A(l,l)=-3/(2*delta_y);
        A(l,l-mnode)=4/(2*delta_y);
        A(l,l-2*mnode)=-1/(2*delta_y);
        r(l)=0;
    end
end
end
end
//Solution of system 
h=A\r
for j=1:nnode
    for i=1:mnode
        l=i+(j-1)*mnode;
        hdata(i,j)=h(l);
    end   
end
contour(x,y,hdata,30)
xtitle("2D model","x axis","y axis");
plot([0 300],[100 100],'-')
plot([300,300],[0 100],'-')
























