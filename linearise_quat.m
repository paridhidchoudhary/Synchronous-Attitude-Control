function [a_quat b_quat]=linearise_quat(A,B,q_0)
E=blkdiag(G(q_0),eye(3));
a_quat=E'*A*E;
b_quat=E'*B;
end