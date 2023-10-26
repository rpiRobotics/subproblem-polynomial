function [x, a, b, c] = half_tan_sp3(p1, p2, k, d)
    [x, a, b, c] = half_tan_sp4(p2, p1, k, 1/2 * (dot(p1,p1)+dot(p2,p2)-d^2));
end