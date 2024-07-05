function [l]=L(q)
s=q(1);
v=q(2:4);
l=[s -v';v s*eye(3,3)*hat(v)];
end