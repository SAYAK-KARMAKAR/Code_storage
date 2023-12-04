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
omega=1.0;
//calculated values of parameters
delta_x=Lx/(mnode-1); 
delta_y=Ly/(nnode-1);
x=0:delta_x:Lx;
y=0:delta_y:Ly;
alphax=1/(delta_x)^2;
alphay=1/(delta_y)^2;
eps_max=1e-3;
//initialization
h=hA*ones(mnode,nnode);
count=0;
rmse=1
while rmse>eps_max
    rmse=0;
    for j=1:nnode
        for i=1:mnode
            if (i>1 & i<mnode) then
                if (j>1 & j<nnode) then //interor nodes
                    cencoff=-2*(alphax+alphay);
                    res=-alphay*h(i,j-1)-alphax*h(i-1,j)+2*(alphax+alphay)*h(i,j)-alphax*h(i+1,j)-alphay*h(i,j+1);
                    h(i,j)=h(i,j)+omega*res/cencoff
                    rmse=rmse+(omega*res/cencoff).^2;
                end
            end
            //Node A
            if (i==1 & j==nnode)then
                h(i,j)=hA;
            end
            //Node B
            if (i==1 & j==1) then h(i,j)=hB; end
            //Node C
            if (i==mnode & j==1) then h(i,j)=hC; end
            //Node D
            if (i==mnode & j==nnode) then h(i,j)=hD;end
            //Specified LBC
            if (i==1) then 
                if (j>1 & j<nnode)then
                    h(i,j)=hB+(hA-hB)*(j-1)*(delta_y/Ly);
                end
            end
            //Specified RBC
            if (i==mnode) then
                if (j>1 & j<nnode) then
                    h(i,j)=hC+(hD-hC)*(j-1)*(delta_y/Ly);
                end
            end
            //Neuman BBC
if (j==1) then
    if(i>1 & i<mnode) then
        h(i,j)=(4*h(i,j+1)-h(i,j+2))/3;
        end
end
//Neumann TBC
if (j==nnode) then
    if (i>1 & i<mnode)then
        //3 point
        h(i,j)=(4*h(i,j-1)-h(i,j-2))/3;
    end
end
end
end
            
rmse=sqrt(rmse/(mnode*nnode)); //overall rmse value
count=count+1
disp([count rmse])
end

contour(x,y,h,50)
xtitle("2D model","x axis","y axis");
plot([0 300],[100 100],'-')
plot([300,300],[0 100],'-')
