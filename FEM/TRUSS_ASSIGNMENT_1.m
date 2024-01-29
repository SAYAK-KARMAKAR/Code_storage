clear all;
clc;
E=200e9*ones(13,1);
A=(1/1000)*[4.5;4.5;4.5;4.5;11.25;7.2;5.625;1.8;5.625;7.2;11.25;9;9];
% area of each element
Node=[0 0;4.8 0;9.6 0;14.4 0;19.2 0;14.4 6.4;9.6 6.4;4.8 6.4];
%coordinates of nodes
elcon=[1 2 E(1) A(1);2 3 E(2) A(2); 3 4 E(3) A(3); 4 5 E(4) A(4);1 8 E(5) A(5);2 8 E(6) A(6); 3 8 E(7) A(7);3 7 E(8) A(8); 3 6 E(9) A(9); 4 6 E(10) A(10);5 6 E(11) A(11);8 7 E(12) A(12); 7 6 E(13) A(13)];

%element connectivity matrix. nodes connected to each element.
UBC=[1 1 0;1 2 0;5 2 0];
% displacement boundary conditions.[node no    x or y   displacement value]
FBC=[2 2 -180; 3 2 -270; 4 2 -360; 5 1 0];
%force boundary condition
numnode=size(Node,1);
%returns the number of rows in the array or matrix.
% Numnode means total number of nodes=8 . if size(Node,2) was there
% it may give number of columns i.e 2
numel=size(elcon,1) ;
%number of elements
Kg=zeros(2*numnode);
%global stiffness matrix is 2 times number of nodes square matrix.
Fg=zeros(2*numnode,1);
%global force vector.
Ug=zeros(2*numnode,1);
%global displacement vector
for el=1:numel
    n1=elcon(el,1);
    n2=elcon(el,2);
    %from elcon matrix, two node numbers are being stored from each element number.
    x1=Node(n1,1);
    y1=Node(n1,2);
    % coordinates of node 1 of an element
    x2=Node(n2,1);
    y2=Node(n2,2);
    % coordinates of node 2 of an element
    theta=atan2(y2-y1,x2-x1);
    %atan2 represents 4 quadrant inverse tangent.
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    C=cos(theta);
    S=sin(theta);
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
Kg([1 2 10],:)=[];
Kg(:,[1 2 10])=[];
K=Kg;
%applying boundary conditions
K;
f=[0;-180000;0;-270000;0;-360000;0;0;0;0;0;0;0];
u=K\f;
displacement=[u(5,1) u(6,1)] %in m

