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
omega=1;
delta_x=Lx/(mnode-1);
delta_y=Ly/(nnode-1);
x=0:delta_x:Lx;
y=0:delta_y:Ly;
delta_t=0.5;//Implicit scheme is unconditionally stable.No need to update
alphax=(T*delta_t)/(S*delta_x^2);
alphay=(T*delta_t)/(S*delta_y^2);
//Initialization
ho=hA*ones(mnode,nnode);
hn=zeros(mnode,nnode);
//Boundary Condition
for j=1:nnode
    //Dirichlet LBC
    ho(1,j)=hB+(hA-hB)*(j-1)*(delta_y/Ly);
    //Dirichlet RBC
    ho(mnode,j)=hC+(hD-hC)*(j-1)*(delta_y/Ly);
end
//Time loop
t=0;
while t<Time_max
    t=t+delta_t;
    count=0;
    rmse=1;
    //space loop
    while rmse>eps_max
        rmse=0;
        for j=1:nnode
           
        for i=1:mnode
            if (i>1 & i<mnode) then
                if (j>1 & j<nnode) then //interor nodes
                    cencoff=-(1+2*(alphax+alphay));
                    res=-ho(i,j)-alphay*hn(i,j-1)-alphax*hn(i-1,j)+(1+2*(alphax+alphay)*hn(i,j))-alphax*hn(i+1,j)-alphay*hn(i,j+1);
                    hn(i,j)=hn(i,j)+omega*res/cencoff
                    rmse=rmse+(omega*res/cencoff).^2;
                end
            end
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
            if(j==1)then
                if (i>1 & i<mnode) then
                    res=hn(i,j+1)-hn(i,j);
                    hn(i,j)=hn(i,j)+omega*(hn(i,j+1)-hn(i,j));//evaluated considering 2 points
                    rmse=rmse+(omega*res).^2; //sum of individual errors
                end
            end
            //Neuman TBC
            if(j==nnode)then
                if (i>1 & i<mnode) then
                    res=hn(i,j-1)-hn(i,j);
                    hn(i,j)=hn(i,j)+omega*(hn(i,j-1)-hn(i,j));
                    rmse=rmse+(omega*res).^2;
                end
            end
        end
    end  
    rmse=sqrt(rmse/(mnode*nnode)); //overall rmse value
count=count+1
disp([count rmse])
end
rmse=0;
for j=1:nnode
    for i=1:mnode
        rmse=rmse+(hn(i,j)-ho(i,j)).^2;
        ho(i,j)=hn(i,j);
end
end

if (rmse<eps_max) then
    break         //while loop terminated
end
end
contour(x,y,hn,30)
xtitle("2D model","x axis","y axis");
plot([0 300],[100 100],'-')
plot([300,300],[0 100],'-')
