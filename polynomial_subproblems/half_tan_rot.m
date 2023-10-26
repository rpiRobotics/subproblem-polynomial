function R = half_tan_rot(k, x)
    s = 2*x/(x^2+1);
    c = (-x^2+1)/(x^2+1);
    R = k*k' + s*hat(k) - c*hat(k)^2;
end