function [x, a, b, c] = half_tan_sp4(h, p, k, d)
    a = 2*h'*k*k'*p - d - h'*p;
    b = 2*h'*cross(k,p);
    c = h'*p-d;
    x = -b +[1 -1]*sqrt(b^2-4*a*c);
    x = x / (2*a);
end