function x = half_tan_sp1(p1, p2, k)
    num = -dot(p2, cross(k,p1));
    den = 2*dot(p1,k)^2 - dot(p1,p2) - dot(p1,p1);
    % den = -dot(p2, cross(k,cross(k,p1)));
    x = num/den;
end