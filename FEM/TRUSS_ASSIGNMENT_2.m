clear all;
clc;
E=10^9*[200; 200; 200;200;200];
A=(1/1000)*[6.4; 6.4; 4; 2.4;4];
% area of each element
Node=[0,0;4.5,3.6;9,0;4.5,6];
%coordinates of nodes
elcon=[1,2;2,3;3 4;2 4;1 4];
%element connectivity matrix. nodes connected to each element.
UBC=[1 1 0;1 2 0; 3 2 0];
% displacement boundary conditions.[node no    x or y   displacement value]
FBC=[4 1 54000];
%force boundary condition
numnode=size(Node,1);
%returns the number of rows in the array or matrix.
% Numnode means total number of nodes=8 . if size(Node,2) was there
% it may give number of columns i.e 2
numel=size(elcon,1); 
%number of elements
Kg=zeros(2*numnode);
%global stiffness matrix is 2 times number of nodes square matrix.
Fg=zeros(2*numnode,1);
%global force vector.
Ug=zeros(2*numnode,1);
%global displacement vector
for el=1:numel
    n1=elcon(el,1);
    n2=elcon(el,2) ;
    %from elcon matrix, two node numbers are being stored from each element number.
    x1=Node(n1,1);
    y1=Node(n1,2);
    % coordinates of node 1 of an element
    x2=Node(n2,1);
    y2=Node(n2,2);
    % coordinates of node 2 of an element
    theta=atan2d(y2-y1,x2-x1);
    %atan2 represents 4 quadrant inverse tangent.
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    C=cosd(theta);
    S=sind(theta);
    kel=(A(el,1)*E(el,1)/L)*[C^2 C*S -C^2 -C*S;C*S S^2 -C*S -S^2;-C^2 -C*S C^2 C*S;-C*S -S^2 C*S S^2];
k1=2*n1-1;
k2=2*n1;
k3=2*n2-1;
k4=2*n2;
%k1, k2,k3,k4 represents the positions of displacement terms of elemental
%stiffness matrix in global stiffness matrix
Kg(k1:k2,k1:k2)=Kg(k1:k2,k1:k2)+kel(1:2,1:2);
Kg(k1:k2,k3:k4)=Kg(k1:k2,k3:k4)+kel(1:2,3:4);
Kg(k3:k4,k1:k2)=Kg(k3:k4,k1:k2)+kel(3:4,1:2);
Kg(k3:k4,k3:k4)=Kg(k3:k4,k3:k4)+kel(3:4,3:4);
end
Kg([1 2 6],:)=[];
Kg(:,[1 2 6])=[];
K=Kg;
%applying boundary conditions
K;
f=[0;0;0;54000;0];
u=K\f;
ux_B=u(3) % Horizontal displacement at B in m

