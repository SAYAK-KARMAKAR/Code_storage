%% Informations given
L=24;
E=200e6;
I=(0.4)^4/12;
F01=36;%Forces are in KN
F02=54;
Node=(0:1.5:L)';% coordinates for nodes
% element connectivity matrix
elcon=zeros(16,2);
for i=1:16
    elcon(i,1)=i;
    elcon(i,2)=i+1;
end
%%
numnode=size(Node,1); %number of nodes
numel=size(elcon,1); %number of elements
Kg=zeros(2*numnode);
Fg=zeros(2*numnode,1);
Ug=zeros(2*numnode,1);
for el=1:numel
    n1=elcon(el,1); %1st node for an element
    n2=elcon(el,2); %2nd node for an element
    le=L/numel;
    kel=((E*I)/le^3)*[12 6*le -12 6*le; 6*le 4*le^2 -6*le 2*le^2; -12 -6*le 12 -6*le; 6*le 2*le^2 -6*le 4*le^2];
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
Kg;
Kg([1 13 33],:)=[];
Kg(:,[1 13 33])=[];
Fg(9,1)=-36;
Fg(25,1)=-54;
Fg([1 13 33],:)=[];
u=Kg\Fg;
U_y=1000*u([6,15,21],:) %Vertical displacements at 4.5m, 12m, 16.5m.
U_theta=(180/pi)*u([1,12,31],:) %Rotation at A,B and C.
%% SFD & BMD using FEM
%Getting total displacement vector
u_total=zeros(2*numnode,1);
u_total(2:12,:)=u(1:11,:);
u_total(14:32)=u(12:30,:);
u_total(34,1)=u(31,1);
%To store SF and BM in "F_EL" for each node
F_EL=[];
for i=1:numel
    u_elemental=[u_total(2*i-1);u_total(2*i);u_total(2*i+1);u_total(2*i+2)];
    f_el=kel*u_elemental;
    F_EL=[F_EL;f_el];
end
for i=1:numel-1 %Because, last two rows for V and M would not come twice.No need to delete.
    F_EL([2*i+1 2*i+2],:)=[];
end
%To extract SF and BM vector seperately from "F_EL"
BM=[];
SF=[];
for i=1:numnode
    BM=[BM;F_EL(2*i)];
    SF=[SF;F_EL(2*i-1)];
end
%Plotting
subplot(2,1,1)
plot(Node,1000*SF);
xlabel('Distance in m');
ylabel('Shear Force in N');
title('SFD');
grid on;
subplot(2,1,2)
plot(Node,BM);
xlabel('Distance in m');
ylabel('Bending moment in kN.m');
title('BMD');
grid on;

%SFD AND BMD Using CLAPEYRON's three moment theorem
% %%  BMD calculations
% % The points where hinge supports are located are named
% %as A,B,C. The points where 36kN and 54 kN loads are acting
% % are named as D and E.
% BMSSB_D=(F01*(Node(5)-Node(1))*(Node(7)-Node(5)))/(Node(7)-Node(1));
%                        % Bending moment at point D for SSB 'AB'
% BMSSB_E=(F02*(Node(13)-Node(7))*(Node(17)-Node(13)))/(Node(17)-Node(7));
%                        % Bending moment at point E for SSB 'BC'
% M_A=0; %Moment at A
% M_C=0; %Moment at C
% L1=(Node(7)-Node(1)); %Lenght of AB
% L2=(Node(17)-Node(7)); %Lenght of BC
% a1=0.5*L1*BMSSB_D; %Area of BMD for simply supported span AB
% x1_bar=(Node(1)+Node(5)+Node(7))/3; %BMD CG distance for the span AB
% a2=0.5*L2*BMSSB_E; %Area of BMD for simply supported span BC
% x2_bar=(Node(7)+Node(13)+Node(17)-3*Node(7))/3; %Because we need distance from
%                                                 %point B
% %Three moment equation
% Term_1=M_A*L1;
% Term_2=M_C*L2;
% Term_3=(6*a1*x1_bar)/L1;
% Term_4=(6*a2*x2_bar)/L2;
% M_B=(Term_3+Term_4-Term_1-Term_2)/(2*(L1+L2)); %Moment at B
% 
% % Initialization
% BM_L1=zeros(numnode,1);
% BM_L2=zeros(numnode,1);
% BM_L=zeros(numnode,1);
% % BM for AB span (Simply supported)
% for i=1:5
%     BM_L1(i,1)=(BMSSB_D/(Node(5)-Node(1)))*(Node(i)-Node(1));
% end
% for i=5:7
%     BM_L1(i,1)=BMSSB_D-(BMSSB_D/(Node(7)-Node(5)))*(Node(i)-Node(5));
% end
% % BM for BC span (Simply supported)
% for i=7:13
%     BM_L2(i,1)=(BMSSB_E/(Node(13)-Node(7)))*(Node(i)-Node(7));
% end
% for i=13:numnode
%     BM_L2(i,1)=BMSSB_E-(BMSSB_E/(Node(numnode)-Node(13)))*(Node(i)-Node(13));
% end
% % BM for M_B
% for i=1:7
%     BM_L(i,1)=(M_B/(Node(7)-Node(1)))*(Node(i)-Node(1));
% end
% for i=7:numnode
%     BM_L(i,1)=M_B-(M_B/(Node(numnode)-Node(7)))*(Node(i)-Node(7));
% end
% BM_L1L2=BM_L1+BM_L2;
% BM_FINAL=BM_L1L2-BM_L;
% 
% % SHEAR FORCE CALCULATIONS
% R_A=(M_B+F01*(Node(7)-Node(5)))/(Node(7)-Node(1));
%                     %Support reaction at A.(Obtained from moment equilibrium
%                     %wrt point B for span AB.
% R_C=(M_B+F02*(Node(13)-Node(7)))/(Node(17)-Node(7));
%                      %Support reaction at C.(Obtained from moment equilibrium
%                     %wrt point B for span BC.
% R_B=(F01+F02)-(R_A+R_C); %Equilibrium of vertical force
% SHEAR_FORCE=(x>=0 & x<6).*(R_A)+(x>=6 & x<9).*(R_A-F01)+(x>=9 & x<18).*(R_A-F01+R_B)+(x>=18 & x<=24).*(R_A-F01+R_B-F02);
% x=linspace(0,24,1000);
% subplot(2,1,1)
% plot(x,SHEAR_FORCE);
% xlabel('Distance in m');
% ylabel('Shear Force in kN');
% title('SFD');
% grid on;
% subplot(2,1,2)
% plot(Node,BM_FINAL);
% xlabel('Distance in m');
% ylabel('Bending moment in kN.m');
% title('BMD');
% grid on;
