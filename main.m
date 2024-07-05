%main
m=4;
w=0.1;
h=0.1;
d=0.1;

b1=0.3;
b2=0.8;

c1=1; %capacitance of battery of cubesats
c2=1;

%moment of inertias of individual cubesats
i1=(m*(h^2 + d^2)/12)
i2=(m*(h^2 + w^2)/12)
i3=(m*(w^2 + d^2)/12)

%moment of inertia post docking
i11=2*m*(h^2 + d^2)/12
i22=2*m*(h^2 + (2*w)^2)/12
i33=2*m*((2*w)^2 + d^2)/12

% i1= 5;
% i2= 6;
% i3= 7;

%equilibrium points
omega_1=0;
omega_2=0;
omega_3=0;
q1=0.1;
q2=0.05;
q3=0.03;
h1=0;
h2=0;
h3=0;

%state space and gains for individual cubesats
% [A,B]=dynamics(omega_1,omega_2,omega_3,i1,i2,i3,q1,q2,q3,q4,h1,h2,h3)
[A,B]=dynamics_new(omega_1,omega_2,omega_3,q1,q2,q3,i1,i2,i3)
% C=[0 0 0 1 0 0;
%     0 0 0 0 1 0;
%     0 0 0 0 0 1];
% D=zeros(3, 3);
C=eye(length(A));
D=zeros(length(A),3);
% [a_quat b_quat]=linearise_quat(A,B,[q1;q2;q3;q4])
is_controlable=check_controllability(A,B);
K=control_sys(A,B);
kd=K(:,1:3);
kp=K(:,4:6);

%overall dynamics
[A_ext,B_ext]= dynamics_ext(i1,i2,i3,i1,i2,i3,i11,i22,i33,q1,q2,q3,q1,q2,q3,q1,q2,q3,omega_1,omega_2,omega_3,omega_1,omega_2,omega_3,omega_1,omega_2,omega_3,c1,c2)
% C_ext = zeros(6, 18);
% C_ext(1:6, 13:18) = eye(6);
% D_ext = zeros(6, 6);
% C_ext=eye(length(A_ext));
% D_ext=zeros(length(A_ext),6);
% 
% Ts=0.2;
% sys_ext=ss(A_ext,B_ext,C_ext,D_ext);
% sysd_ext=c2d(sys_ext,Ts);
% mpcobj= mpc(sysd_ext);
% 
% mpcobj.PredictionHorizon=10;
% mpcobj.ControlHorizon=3;
% 
% for i=1:6
%     mpcobj.ManipulatedVariables(i).Min=-0.25;
%     mpcobj.ManipulatedVariables(i).Max=0.25;
% end
% 
% for i=[1:3,7:9,13:15]
%     mpcobj.OutputVariables(i).Min= -0.1;
%     mpcobj.OutputVariables(i).Max= 0.1;
% end



