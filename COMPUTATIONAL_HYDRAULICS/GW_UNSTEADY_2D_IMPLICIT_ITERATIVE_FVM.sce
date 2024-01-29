clc
clear
//GW unsteady 2D iterative implicit FVM
//....Problem dependent parameters..
mnode=31;
nnode=21;
mcell=mnode-1;
ncell=nnode-1;
Lx=300;
Ly=100;
hA=90;
hB=89;
hC=85;
hD=87;
Time_max=50;
eps_max=1e-3;
S=5e-5;//storativity
T=200;//Transmissivity
omega=1;
delta_x=Lx/(mnode-1);
delta_y=Ly/(nnode-1);
x=delta_x/2:delta_x:(Lx-delta_x/2);
y=delta_y/2:delta_y:(Ly-delta_y/2);
delta_t=0.5;
alphax=(T*delta_t)/(S*delta_x^2);
alphay=(T*delta_t)/(S*delta_y^2);
//modified corner heads (Because, these values should be defined at cell centered positions, not at the corners anymore)
hAm=hB+((hA-hB)/Ly)*(Ly-delta_y/2); //Ly-delta_y/2 is the distance between corner point B and cell centered A.
hBm=hB+((hA-hB)/Ly)*(delta_y/2);
hCm=hC+((hD-hC)/Ly)*(delta_y/2);
hDm=hC+((hD-hC)/Ly)*(Ly-delta_y/2);
//Initialization 
ho=hA*ones(mcell,ncell);//old time level
hn=hA*ones(mcell,ncell);//New time level
//Time loop
t=0;
while t<Time_max
    t=t+delta_t;
    count=0;
    rmse=1;
//space loop
while rmse>eps_max
    rmse=0;
    for j=1:ncell
        for i=1:mcell
//Interior nodes
if(i>1 & i<mcell) then
    if(j>1 & j<ncell)then
        a_S=alphay;
        a_W=alphax;
        a_P=-1-(2*(alphax+alphay));
        a_E=alphax;
        a_N=alphay;
        r_P=-ho(i,j);
        res=r_P-(a_S*hn(i,j-1)+a_W*hn(i-1,j)+a_P*hn(i,j)+a_E*hn(i+1,j)+a_N*hn(i,j+1));
end
end
//cell A
if (i==1 & j==ncell) then
   a_S=alphay;
   a_W=0;
   a_P=-1-4*alphax-alphay;
   a_E=(4/3)*alphax;
   a_N=0; 
end
//Cell B
if (i==1 & j==1) then
    a_S=0;
    a_W=0;
    a_P=-1-4*alphax-alphay;
    a_E=(4/3)*alphax;
    a_N=alphay;
    r_P=-ho(i,j)-(8/3)*alphax*hBm;
    res=r_P-(a_P*hn(i,j)+a_E*hn(i+1,j)+a_N*hn(i,j+1));
end
//Cell C
if (i==mcell & j==1) then
    a_P=-1-4*alphax-alphay;
    res=(-ho(i,j)-(8/3)*alphax*hCm)-((4/3)*alphax*hn(i-1,j)-(1+4*alphax+alphay)*hn(i,j)+alphay*hn(i,j+1));
end
//Cell D
if (i==mcell & j==ncell) then
    a_P=-(1+4*alphax+alphay);
    res=(-ho(i,j)-(8/3)*alphax*hDm)-(alphay*hn(i,j-1)+(4/3)*alphax*hn(i-1,j)-(1+4*alphax+alphay)*hn(i,j));
end
//Specified LBC
if (i==1) then
    if (j>1 & j<ncell)then
      a_P=-(1+4*alphax+2*alphay);  
    hvl=hBm+(hAm-hBm)*(j-1)*((delta_y)/(Ly-delta_y)); //hvl=head value at the left boundary.As FVM considers cell centered values dimensions of the domain length=(Lx-delta_x) and width=(Ly-delta_y)
    res=(-ho(i,j)-(8/3)*alphax*hvl)-(alphay*hn(i,j-1)-(1+4*alphax+2*alphay)*hn(i,j)+(4/3)*alphax*hn(i+1,j)+alphay*hn(i,j+1));
end
end
//Specified RBC
if (i==mcell) then
    if (j>1 & j<ncell)then
        a_P=-(1+4*alphax+2*alphay);
        hvr=hCm+(hDm-hCm)*(j-1)*delta_y/(Ly-delta_y);
        res=(-ho(i,j)-(8/3)*alphax*hvr)-(alphay*hn(i,j-1)-(1+4*alphax+2*alphay)*hn(i,j)+(4/3)*alphax*hn(i-1,j)+alphay*hn(i,j+1))
end
//Neumann TBC
if (j==ncell) then
    if (i>1 & i<mcell)then
        a_P=-(1+2*alphax+alphay);
        res=-ho(i,j)-(alphay*hn(i,j-1)+alphax*hn(i-1,j)-(1+2*alphax+alphay)*hn(i,j)+alphax*hn(i+1,j));
end
end
//Neumann BBC
if (j==1) then
    if (i>1 & i<mcell)then
        a_P=-(1+2*alphax+alphay);
        res=-ho(i,j)-(alphax*hn(i-1,j)-(1+2*alphax+alphay)*hn(i,j)+alphax*hn(i+1,j)++alphay*hn(i,j+1));
end
end
//Update
hn(i,j)=hn(i,j)+omega*res/a_P;
rmse=rmse+(omega*res/a_P).^2;//Total rmse
end
end
rmse=sqrt(rmse/(mcell*ncell));//Actual rmse
count =count +1;
disp([count rmse])
end
//To check whether things are converging with time or not
rmse_t=0;
for j=1:ncell 
    for i=1:mcell
        rmse_t=rmse_t+(hn(i,j)-ho(i,j)).^2;
        ho(i,j)=hn(i,j);
    end
end
//Condition for steady state
if (rmse_t<eps_max) then
    break
end
end
//Boundary information
hdata=zeros(ncell+2,mcell+2);//corner points need to be added to this matrix to get the picture of full domain. those points were excluded during FVM formulation.
//Internal cells
for j=2:ncell+1;
    for i=mcell+1
        hdata(i,j)=hn(i-1,j-1);//because, after updating, one extra row of nodes added to the bottom.
    end
end
//Input the corner values at 'hdata' matrix
hdata(1,ncell+2)=hA;
hdata(1,1)=hB;
hdata(mcell+2,1)=hC;
hdata(mcell+2,ncell+2)=hD;
//Left and right boundary
for j=2:ncell+1
    hdata(1,j)=hBm+(hAm-hBm)*((j-2)*delta_y)/(Ly-delta_y);//left boundary. along a column of 'hdata',1st node,h(1,1)=hB; 2nd node h(2,1)=hBm; etc. 

hdata(mcell+2,j)=hCm+(hDm-hCm)*(j-2)*(delta_y/(Ly-delta_y)); //Right boundary 
end
//Bottom and top boundary
for i=2:mcell+1
    hdata(i,1)=(9*hn(i-1,1)-hn(i-1,2))/8;//Bottom boundary
    hdata(i,ncell+2)=(9*hn(i-1,1)-hn(i-1,ncell-1))/8;
    
end
end
xn=[0 x Lx];
yn=[0 y Ly];
contour(xn,yn,hdata,30)
xtitle("2D FVM","x axis","y axis");
plot([0 300],[100 100],'-')
plot([300,300],[0 100],'-')
    

















