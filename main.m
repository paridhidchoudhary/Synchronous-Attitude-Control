%main
m=4;
w=0.3;
h=0.1;
d=0.1;

b1=0.3;
b2=0.8;

% i1=(m*(h^2 + d^2)/12)
% i2=(m*(h^2 + w^2)/12)
% i3=(m*(w^2 + d^2)/12)

i1= 5;
i2= 6;
i3= 7;

%equilibrium points
omega_1=0.1;
omega_2=0.1;
omega_3=0.1;
q1=1;
q2=0;
q3=0;
q4=0;
h1=0;
h2=0;
h3=0;

[A,B]=dynamics(omega_1,omega_2,omega_3,i1,i2,i3,q1,q2,q3,q4,h1,h2,h3)
C=eye(length(A));
D=zeros(7,3);
is_controlable=check_controllability(A, B)
K=control_sys(A,B)
kp=[10 0 0 0;0 0.0001 0 0;0 0 0 0.0001]
kd=[0.01 0 0;0 0.01 0;0 0 0.01]