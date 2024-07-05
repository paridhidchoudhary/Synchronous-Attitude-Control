function [A_sys,B_sys]= dynamics(omega_1_eq,omega_2_eq,omega_3_eq,I1_sys,I2_sys,I3_sys,q1_eq,q2_eq,q3_eq,q4_eq,h1_eq,h2_eq,h3_eq)
syms omega_1 omega_2 omega_3 I1 I2 I3 tau_1 tau_2 tau_3 q1 q2 q3 q4 h1 h2 h3;
% q=[q1; q2; q3; q4];
g=[q1;q2;q3];
w=[omega_1;omega_2;omega_3];
%h=[h1;h2;h3];
omega=[0 omega_3 -omega_2 omega_1;-omega_3 0 omega_1 omega_2;omega_2 -omega_1 0 omega_3;-omega_1 -omega_2 -omega_3 0];
I=[I1 0 0;0 I2 0; 0 0 I3];

tau=[tau_1;tau_2;tau_3];
omega_skew=[0 -omega_3 omega_2;omega_3 0 -omega_1;-omega_2 omega_1 0];

% dqdt=0.5 * omega* q;
%dgdt=-0.5*omega_skew*g + 0.5*q4*w;
dgdt=0.5*eye(3)*w;
% dq4dt=-0.5*w'*g;
omega_skew*I*w
I_inv=[1/I1 0 0;0 1/I2 0;0 0 1/I3];
domegadt= I_inv*(tau- (omega_skew*I*w))%-I_inv*omega_skew*h;
%dhdt=-tau;

%f=[dqdt;domegadt;dhdt]
% f=[dqdt; domegadt]
f=[dgdt ;domegadt];
A=jacobian(f, [g;w])
B=jacobian(f, tau)

A_sys=double(subs(A,{omega_1,omega_2,omega_3,I1,I2,I3,q1,q2,q3,q4},{omega_1_eq,omega_2_eq,omega_3_eq,I1_sys,I2_sys,I3_sys,q1_eq,q2_eq,q3_eq,q4_eq}));
B_sys=double(subs(B,{I1,I2,I3},{I1_sys,I2_sys,I3_sys}));

end
