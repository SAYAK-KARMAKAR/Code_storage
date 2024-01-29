clc
clear
//....Problem dependent paramaters...
mnode=31;
nnode=21;
Lx=300;
Ly=100;
hA=90;
hB=89;
hC=85;
hD=87;
Time_max=5;
eps_max=1e-3;
S=5e-5;//storativity
T=200;//Transmissivity
//Calculated parameter values
delta_x=Lx/(mnode-1);
delta_y=Ly/(nnode-1);
x=0:delta_x:Lx;
y=0:delta_y:Ly;
delta_t=0.5;//we assumed. May need to be recalculated based on alphax and alphay values.
alphax=(T*delta_t)/(S*delta_x^2);
alphay=(T*delta_t)/(S*delta_y^2);
sumalpha=alphax+alphay;
while sumalpha>0.5
    delta_t=delta_t/2
    alphax=(T*delta_t)/(S*delta_x^2);
    alphay=(T*delta_t)/(S*delta_y^2);
    sumalpha=alphax+alphay;//this while loop is to ensure stability criteria of explicit problem and to take appropriate delta_t value. 
end
disp("delta_t",delta_t)
//Initialization
ho=hA*ones(mnode,nnode);//ho for old time level
hn=zeros(mnode,nnode);//hn for new time level
//Boundary Condition
for j=1:nnode
    //Dirichlet LBC
    ho(1,j)=hB+(hA-hB)*(j-1)*(delta_y/Ly);
    //Dirichlet RBC
    ho(mnode,j)=hC+(hD-hC)*(j-1)*(delta_y/Ly);
end
count=0;
rmse=1; 
t=0;
//Time loop
while t<Time_max
    t=t+delta_t;
    //Interior node
    for j=1:nnode
        for i=1:mnode
            if (i>1 & i<mnode)then
                if (j>1 & j<nnode)then
                    hn(i,j)=alphay*ho(i,j-1)+alphax*ho(i-1,j)+(1-2*(alphax+alphay))*ho(i,j)+alphax*ho(i+1,j)+alphay*ho(i,j+1);                  
end
end
end
end
//Boundary nodes
for j=1:nnode
    for i=1:mnode
        //node A
            if (i==1 & j==nnode)then
                hn(i,j)=hA;
            end
            //Node B
            if (i==1 & j==1) then hn(i,j)=hB; end
            //Node C
            if (i==mnode & j==1) then hn(i,j)=hC; end
            //Node D
            if (i==mnode & j==nnode) then hn(i,j)=hD;end
            //Specified LBC
            if (i==1) then 
                if (j>1 & j<nnode)then
                    hn(i,j)=hB+(hA-hB)*(j-1)*(delta_y/Ly);
                end
            end
            //Specified RBC
            if (i==mnode) then
                if (j>1 & j<nnode) then
                    hn(i,j)=hC+(hD-hC)*(j-1)*(delta_y/Ly);
                end
            end
//Neuman BBC
if (j==1) then
    if(i>1 & i<mnode) then
        hn(i,j)=(4*hn(i,j+1)-hn(i,j+2))/3;
        end
end
//Neumann TBC
if (j==nnode) then
    if (i>1 & i<mnode)then
        //3 point
        hn(i,j)=(4*hn(i,j-1)-hn(i,j-2))/3; //Boundary conditions are evaluated based on new time level values
end
end
end
end
rmse=0;
for j=1:nnode
    for i=1:mnode
        rmse=rmse+(hn(i,j)-ho(i,j)).^2;
        ho(i,j)=hn(i,j);
end
end
rmse=sqrt(rmse/(mnode*nnode));
disp("t rmse",[t rmse])
if (rmse<eps_max) then
    break         //while loop terminated
end
end
contour(x,y,hn,50)
xtitle("2D model","x axis","y axis");
plot([0 300],[100 100],'-')
plot([300,300],[0 100],'-')






















