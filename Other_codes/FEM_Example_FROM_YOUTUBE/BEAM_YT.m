clear;
clc;
E=200e9;
I=4/10^6
%L=24;
q=[-20000; -10000];
xi=[-sqrt(1/3); sqrt(1/3)];
w=[1; 1];
%% element connectivity matrix
elcon=zeros(2,2);
for i=1:2
    elcon(i,1)=i;
    elcon(i,2)=i+1;
end
numel=size(elcon,1);
le=[4;5]
elcon
numnode=numel+1;
%% coordinates for nodes
Node_coordinate=[0 0;4 0; 9 0];

%% Initialization 
Kg=zeros(2*numnode)
%Ug=zeros(2*numnode,1);
Fg=zeros(2*numnode,1) 
%N1=0.25*(2-3*xi+xi^3)
%N2=(le/8)*(1-xi-xi^2+xi^3)
%N3=0.25*(2+3*xi-xi^3)
%N4=(le/8)*(-1-xi+xi^2+xi^3)
%% for 1st gauss point
ddN1_1=1.5*xi(1);
ddN2_1=(le/8)*(-2+6*xi(1));
ddN3_1=-1.5*xi(1);
ddN4_1=(le/8)*(2+6*xi(1));
ddN_VEC_1=[ddN1_1; ddN2_1; ddN3_1; ddN4_1];
%% for 2nd gauss point
ddN1_2=1.5*xi(2);
ddN2_2=(le/8)*(-2+6*xi(2));
ddN3_2=-1.5*xi(2);
ddN4_2=(le/8)*(2+6*xi(2));
ddN_VEC_2=[ddN1_2; ddN2_2; ddN3_2; ddN4_2];
%% elemental stiffness matrix
kel_1=zeros(4);
kel_2=zeros(4);
for i=1:4
    for j=1:4
kel_1(i,j)=E*I*(le(1)/2)^3*(ddN_VEC_1(i)*ddN_VEC_1(j)*w(1)+ddN_VEC_2(i)*ddN_VEC_2(j)*w(2));
kel_2(i,j)=E*I*(le(2)/2)^3*(ddN_VEC_1(i)*ddN_VEC_1(j)*w(1)+ddN_VEC_2(i)*ddN_VEC_2(j)*w(2));
    end
end
kel_1
kel_2
% %% for 1st gauss point
% N1_1=0.25*(2-3*xi(1)+xi(1)^3);
% N2_1=(le/8)*(1-xi(1)-xi(1)^2+xi(1)^3);
% N3_1=0.25*(2+3*xi(1)-xi(1)^3);
% N4_1=(le/8)*(-1-xi(1)+xi(1)^2+xi(1)^3);
% N_VEC_1=[N1_1;N2_1;N3_1;N4_1];
% %% for 2nd gauss point
% N1_2=0.25*(2-3*xi(2)+xi(2)^3);
% N2_2=(le/8)*(1-xi(2)-xi(2)^2+xi(2)^3);
% N3_2=0.25*(2+3*xi(2)-xi(2)^3);
% N4_2=(le/8)*(-1-xi(2)+xi(2)^2+xi(2)^3);
% N_VEC_2=[N1_2;N2_2;N3_2;N4_2];
% %% elemental force vector for distributed load
% fel=zeros(4,1);
% for i=1:4
% fel(i)=(le/2)*q*[N_VEC_1(i)*w(1)+N_VEC_2(i)*w(2)];
% end
% for el=1:numel % el=element number, numel=total no of elements
%     n1=elcon(el,1);
%     n2=elcon(el,2);
% k1=2*n1-1;
% k2=2*n1;
% k3=2*n2-1;
% k4=2*n2;
% %k1, k2,k3,k4 represents the positions of displacement terms of elemental
% %stiffness matrix in global stiffness matrix
% Kg(k1:k4,k1:k4)=Kg(k1:k4,k1:k4)+kel(1:4,1:4);
% Fg(k1:k4,1)=Fg(k1:k4,1)+fel(1:4,1);
% end
% Fg(5,1)=-36000;
% Fg(13,1)=-54000;
% Kg([1,13,33],:)=[];
% Kg(:,[1,13,33])=[];
% Fg([1,13,33],:)=[];
% Kg;
% Fg;
% Ug=Kg\Fg