function [M Meq] = constraint_1(G,a,y,e)
M = norm(G * a - y) - e;
Meq = [];
end
