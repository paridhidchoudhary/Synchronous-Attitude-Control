function [g]=G(q)
H=[zeros(1,3); eye(3,3)];
g=L(q)*H;
end