function [x1, x2, a1, b1, c1, a2, b2, c2] = half_tan_sp2(p1, p2, k1, k2)
    [x1, a1, b1, c1] = half_tan_sp4(k2, p1, k1, dot(k2,p2));
    [x2, a2, b2, c2] = half_tan_sp4(k1, p2, k2, dot(k1,p1));
end