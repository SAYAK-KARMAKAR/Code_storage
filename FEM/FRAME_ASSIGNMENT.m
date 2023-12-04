clc
clear all
E=200e9;
A=0.09;
Ic=6.75*10^(-4);
I_VEC=[Ic; Ic; 2*Ic; 2*Ic; Ic];% Area moment of inertia
node=[0 0;0 4.5;0,7.5;3,7.5;6,7.5;6,2.5]; % node coordinates
elcon=[1 2;2 3;3 4;4,5; 5,6]; % element connectivity matrix
numel=size(elcon,1); %number of elements
numnode=size(node,1); %number of nodes
k_global=zeros(3*numnode);
for el=1:numel
    n1=elcon(el,1);
    n2=elcon(el,2);
    % n1,n2 represent node numbers for 'el'th element.
    x1=node(n1,1);
    y1=node(n1,2);
    % (x1,y1) coordinates of node 1 of an element
    x2=node(n2,1);
    y2=node(n2,2);
    % coordinates of node 2 of an element
    theta=atan2d(y2-y1,x2-x1);
    %atan2 represents 4 quadrant inverse tangent.
    L=sqrt((x2-x1)^2+(y2-y1)^2); % LENGTH OF THE ELEMENT
    C=cosd(theta);
    S=sind(theta);
    I=I_VEC(el,1); % area moment of inertia for 'el'th element
    T1=(A*E)/L;
    T2=(12*E*I)/L^3;
    T3=(4*E*I)/L;
    T4=(6*E*I)/L^2;
    T5=(2*E*I)/L;
    k_el_lcs=[T1 0 0 -T1 0 0;
        0 T2 T4 0 -T2 T4;
        0 T4 T3 0 -T4 T5;
        -T1 0 0 T1 0 0;
        0 -T2 -T4 0 T2 -T4;
        0 T4 T5 0 -T4 T3]; %elemental stiffness matrix in local
    % co ordinate system
    TRANSFORMATION_MAT=zeros(6);
    TRANSFORMATION_MAT(1:2,1:2)=[C S;-S C];
    TRANSFORMATION_MAT(4:5,4:5)=[C S;-S C];
    TRANSFORMATION_MAT(3,3)=1;
    TRANSFORMATION_MAT(6,6)=1;
    k_el_gcs=TRANSFORMATION_MAT.'*k_el_lcs*TRANSFORMATION_MAT;
    k1=3*n1-2;
    k2=3*n1-1;
    k3=3*n1;
    k4=3*n2-2;
    k5=3*n2-1;
    k6=3*n2;
    k_global(k1:k6,k1:k6)=k_global(k1:k6,k1:k6)+k_el_gcs;
    %Kg(k1:k2,k3:k4)=Kg(k1:k2,k3:k4)+kel(1:2,3:4)
    %Kg(k3:k4,k1:k2)=Kg(k3:k4,k1:k2)+kel(3:4,1:2)
    %Kg(k3:k4,k3:k4)=Kg(k3:k4,k3:k4)+kel(3:4,3:4)
end
f_global=[0; 0; 0; 48000; 0; 0; 0; 0; 0; 0; -96000; 0; 0; 0; 0; 0; 0; 0];
k_global([1 2 17],:)=[];
k_global(:,[1 2 17])=[];
f_global([1 2 17],:)=[];
u=k_global\f_global
Theta_ABCD=[u(1);u(7);u(13);u(15)] %in radian
ux_D=u(14) %in m