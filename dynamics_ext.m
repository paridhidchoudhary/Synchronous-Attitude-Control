function [A_ext,B_ext]=dynamics_ext(i11,i12,i13,i21,i22,i23,i31,i32,i33,g11,g12,g13,g21,g22,g23,g31,g32,g33,w11,w12,w13,w21,w22,w23,w31,w32,w33,c1,c2)
syms I11 I12 I13 I21 I22 I23 I31 I32 I33 q10 q11 q12 q13 q20 q21 q22 q23 q30 q31 q32 q33 omega11 omega12 omega13 omega21 omega22 omega23 omega31 omega32 omega33 u11 u12 u13 u21 u22 u23 b1 b2
Q_1=[q10 -q13 q12;
    q13 q10 -q11;
    -q12 q11 q10];
Q_2=[q20 -q23 q22;
    q23 q20 -q21;
    -q22 q21 q20];
Q_3=[q30 -q33 q32;
    q33 q30 -q31;
    -q32 q31 q30];
u_1=[u11;u12;u13];
u_2=[u21;u22;u23];
u_3=u_1+u_2;
w1=[omega11;omega12;omega13];
w2=[omega21;omega22;omega23];
w3=[omega31;omega32;omega33];
g1=[q11;q12;q13];
g2=[q21;q22;q23];
g3=[q31;q32;q33];

q10_eq=sqrt(1-g11^2 - g12^2 - g13^2);
q20_eq=sqrt(1-g21^2 - g22^2 - g23^2);
q30_eq=sqrt(1-g31^2 - g32^2 - g33^2);
I1=[I11 0 0;
    0 I12 0;
    0 0 I13];
I1_inv=[1/I11 0 0;
    0 1/I12 0;
    0 0 1/I13];
I2=[I21 0 0;
    0 I22 0;
    0 0 I23];
I2_inv=[1/I21 0 0;
    0 1/I22 0;
    0 0 1/I23];
I3=[I31 0 0;
    0 I32 0;
    0 0 I33];
I3_inv=[1/I31 0 0;
    0 1/I32 0;
    0 0 1/I33];

p_1= u11*omega11 + u12*omega12 + u13*omega13;
p_2= u21*omega21 + u22*omega22 + u23*omega23;

domegadt1=I1_inv*(-hat(w1)*I1*w1 -u_1);
dgdt1=0.5*Q_1*w1;
domegadt2=I2_inv*(-hat(w2)*I2*w2 -u_2);
dgdt2=0.5*Q_2*w2;
domegadt3=I3_inv*(-hat(w3)*I3*w3 -u_3);
dgdt3=0.5*Q_3*w3;
dbdt1= -p_1/c1;
dbdt2= -p_2/c2;

f=[domegadt1;dgdt1;domegadt2;dgdt2;domegadt3;dgdt3;dbdt1;dbdt2];

A=jacobian(f,[w1;g1;w2;g2;w3;g3;b1;b2]);
B=jacobian(f,[u_1;u_2]);

A_ext=double(subs(A,{omega11,omega12,omega13,omega21,omega22,omega23,omega31,omega32,omega33,q10,q11,q12,q13,q20,q21,q22,q23,q30,q31,q32,q33,I11,I12,I13,I21,I22,I23,I31,I32,I33,u11,u12,u13,u21,u22,u23},{w11,w12,w13,w21,w22,w23,w31,w32,w33,q10_eq,g11,g12,g13,q20_eq,g21,g22,g23,q30_eq,g31,g32,g33,i11,i12,i13,i21,i22,i23,i31,i32,i33,0,0,0,0,0,0}));
B_ext=double(subs(B,{I11,I12,I13,I21,I22,I23,I31,I32,I33,u11,u12,u13,u21,u22,u23,omega11,omega12,omega13,omega21,omega22,omega23,omega31,omega32,omega33},{i11,i12,i13,i21,i22,i23,i31,i32,i33,0,0,0,0,0,0,w11,w12,w13,w21,w22,w23,w31,w32,w33}));

end