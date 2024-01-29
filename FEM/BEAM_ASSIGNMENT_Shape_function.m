clear all;
clc;
E=200e9;
I=0.4^4/12;
L=24;
q=0;
xi=[-sqrt(1/3); sqrt(1/3)];
w=[1; 1];
%% element connectivity matrix
elcon=zeros(16,2);
for i=1:16
    elcon(i,1)=i;
    elcon(i,2)=i+1;
end
numel=size(elcon,1);
le=L/numel;
elcon;
numnode=numel+1;
%% coordinates for nodes
Node_coordinate=zeros(numnode,1);
for i=1:numnode
    Node_coordinate(i,1)=Node_coordinate(i,1)+(i-1)*le;
end
Node_coordinate;
%% Initialization 
Kg=zeros(2*numnode);
%Ug=zeros(2*numnode,1);
Fg=zeros(2*numnode,1);
%% Shape Functions
%N1=0.25*(2-3*xi+xi^3)
%N2=(le/8)*(1-xi-xi^2+xi^3)
%N3=0.25*(2+3*xi-xi^3)
%N4=(le/8)*(-1-xi+xi^2+xi^3)
%% 2nd derivative of shape functions for 1st gauss point
ddN1_1=1.5*xi(1);
ddN2_1=(le/8)*(-2+6*xi(1));
ddN3_1=-1.5*xi(1);
ddN4_1=(le/8)*(2+6*xi(1));
ddN_VEC_1=[ddN1_1; ddN2_1; ddN3_1; ddN4_1];
%% 2nd derivative of shape functions for 2nd gauss point
ddN1_2=1.5*xi(2);
ddN2_2=(le/8)*(-2+6*xi(2));
ddN3_2=-1.5*xi(2);
ddN4_2=(le/8)*(2+6*xi(2));
ddN_VEC_2=[ddN1_2; ddN2_2; ddN3_2; ddN4_2];
%% elemental stiffness matrix
kel=zeros(4);
for i=1:4
    for j=1:4
kel(i,j)=E*I*(2/le)^3*(ddN_VEC_1(i)*ddN_VEC_1(j)*w(1)+ddN_VEC_2(i)*ddN_VEC_2(j)*w(2));
    end
end
kel;
%% for 1st gauss point
N1_1=0.25*(2-3*xi(1)+xi(1)^3);
N2_1=(le/8)*(1-xi(1)-xi(1)^2+xi(1)^3);
N3_1=0.25*(2+3*xi(1)-xi(1)^3);
N4_1=(le/8)*(-1-xi(1)+xi(1)^2+xi(1)^3);
N_VEC_1=[N1_1;N2_1;N3_1;N4_1];
%% for 2nd gauss point
N1_2=0.25*(2-3*xi(2)+xi(2)^3);
N2_2=(le/8)*(1-xi(2)-xi(2)^2+xi(2)^3);
N3_2=0.25*(2+3*xi(2)-xi(2)^3);
N4_2=(le/8)*(-1-xi(2)+xi(2)^2+xi(2)^3);
N_VEC_2=[N1_2;N2_2;N3_2;N4_2];
%% elemental force vector for distributed load
fel=zeros(4,1);
for i=1:4
fel(i)=(le/2)*q*(N_VEC_1(i)*w(1)+N_VEC_2(i)*w(2));
end
for el=1:numel % el=element number, numel=total no of elements
    n1=elcon(el,1);
    n2=elcon(el,2);
k1=2*n1-1;
k2=2*n1;
k3=2*n2-1;
k4=2*n2;
%k1, k2,k3,k4 represents the positions of displacement terms of elemental
%stiffness matrix in global stiffness matrix
Kg(k1:k4,k1:k4)=Kg(k1:k4,k1:k4)+kel(1:4,1:4);
Fg(k1:k4,1)=Fg(k1:k4,1)+fel(1:4,1);
end
Fg(5,1)=-36000;
Fg(13,1)=-54000;
Kg([1,13,33],:)=[];
Kg(:,[1,13,33])=[];
Fg([1,13,33],:)=[];
Kg;
Fg;
Ug=Kg\Fg;
U_y=Ug([6,15,21],:) %Vertical displacements at 4.5m, 12m, 16.5m in m.
U_theta=(180/pi)*Ug([1,12,31],:) %Rotation at A,B and C IN DEGREES.



















% %%
% L = 24; % length of beam in meters
% n = 3; % number of supports
% support_pos = [0 9 24]; % positions of the supports in meters
% w1 = 36000; % load in N
% w2 = 54000; % load in N
% load_pos = [6 18]; % positions of the loads in meters
% R = zeros(n+1,1);
% M = zeros(n+1,1);
% R(1) = (w1 + w2)*L/2/(support_pos(2)-support_pos(1)); % reaction at first support
% R(end) = (w1 + w2)*L/2/(support_pos(3)-support_pos(2)); % reaction at last support
% M(1) = 0;
% M(end) = 0;
% for i = 2:n
%     R(i) = R(i-1) - (w1 + w2)*L/2/(support_pos(i)-support_pos(i-1));
%     M(i) = M(i-1) + R(i-1)*(support_pos(i)-support_pos(i-1));
% end
% x = linspace(0,L,241);
% V = zeros(size(x));
% M = zeros(size(x));
% for i = 1:length(x)
%     if x(i) < support_pos(2)
%         V(i) = R(1) - w1*(x(i)-0) - w2*(x(i)-0);
%         M(i) = R(1)*x(i) - w1*(x(i)-0)^2/2 - w2*(x(i)-0)^2/2;
%     elseif x(i) < support_pos(3)
%         V(i) = R(2) - w1*(x(i)-support_pos(2)) - w2*(x(i)-0);
%         M(i) = R(2)*(x(i)-support_pos(2)) + M(2) - w1*(x(i)-support_pos(2))^2/2 - w2*(x(i)-0)^2/2;
%      else
%         V(i) = R(3) - w1*(x(i)-support_pos(2)) - w2*(x(i)-support_pos(3));
%         M(i) = R(3)*(x(i)-support_pos(3)) + M(3) - w1*(x(i)-support_pos(2))^2/2 - w2*(x(i)-support_pos(3))^2/2;
%     end
% end
% figure;
% subplot(2,1,1); % subplot for shear force
% plot(x,V);
% title('Shear Force Diagram');
% xlabel('Distance (m)');
% ylabel('Shear Force (N)');
% grid on;
% subplot(2,1,2); % subplot for bending moment
% plot(x,M);
% title('Bending Moment Diagram');
% SF_VAL=[V(1); V(91); V(241); V(61); V(181)]