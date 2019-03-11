function [q, fjac] = simpleFunc(p)

q = -12 + 2*p.^(-3);

fjac = -6*p.^(-4);

end