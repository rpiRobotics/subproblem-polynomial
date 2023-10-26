function [lhs_1, lhs_2, R] = half_tan_sp5_R(p0, p1, p2, p3, k1, k2, k3, x_vec)
    % sp_6(H, K, P, d1, d2)
    
    H = [p0 -p2 k2 -k2];
    K = [k1 k3 k1 k3];
    P = [p1 p3 p1 p3];
    d1 = 0.5*(dot(p2,p2) + dot(p3,p3) - dot(p0,p0) - dot(p1,p1));
    d2 = dot(k2,p2) - dot(k2,p0);
    [lhs_1, lhs_2] = half_tan_sp6(H, K, P, d1, d2, [x_vec(1) x_vec(2)]);

    R_1 = half_tan_rot(k1, x_vec(1));
    R_3 = half_tan_rot(k3, x_vec(2));
    R = half_tan_sp1_R(p2+R_3*p3, p0+R_1*p1, k2);
    %x2 = half_tan_sp1(p2+R_3*p3, p0+R_1*p1, k2);
    %R = half_tan_rot(k2, x2);
end