function[A_sys, B_sys]=dynamics_new(omega1_eq,omega2_eq,omega3_eq,q1_eq,q2_eq,q3_eq,I1_sys,I2_sys,I3_sys)
syms I1 I2 I3 q0 q1 q2 q3 omega1 omega2 omega3 u1 u2 u3
Q=[q0 -q3 q2;
    q3 q0 -q1;
    -q2 q1 q0];
u=[u1;u2;u3];
w=[omega1;omega2;omega3];
g=[q1;q2;q3];
% q0=sqrt(1-q1^2-q2^2-q3^2);
q0_eq=sqrt(1-q1_eq^2-q2_eq^2-q3_eq^2)
I=[I1 0 0;
    0 I2 0;
    0 0 I3];
I_inv=[1/I1 0 0;
    0 1/I2 0;
    0 0 1/I3];
omega_skew=[0 -omega3 omega2;
    omega3 0 -omega1;
    -omega2 omega1 0];

domegadt=I_inv*(-omega_skew*I*w -u);
dgdt=0.5*Q*w;
f=[domegadt;dgdt];

A=jacobian(f,[w;g]);
B=jacobian(f,u);

A_sys=double(subs(A,{omega1,omega2,omega3,q0,q1,q2,q3,I1,I2,I3},{omega1_eq,omega2_eq,omega3_eq,q0_eq,q1_eq,q2_eq,q3_eq,I1_sys,I2_sys,I3_sys}));
B_sys=double(subs(B,{I1,I2,I3},{I1_sys,I2_sys,I3_sys}));

end
