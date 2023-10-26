function [lhs_1, lhs_2] = half_tan_sp6(H, K, P, d1, d2, x_vec)
R_1 = half_tan_rot(K(:,1), x_vec(1));
R_2 = half_tan_rot(K(:,2), x_vec(2));
R_3 = half_tan_rot(K(:,3), x_vec(1));
R_4 = half_tan_rot(K(:,4), x_vec(2));

lhs_1 = H(:,1)'*R_1*P(:,1) + H(:,2)'*R_2*P(:,2) - d1;
lhs_2 = H(:,3)'*R_3*P(:,3) + H(:,4)'*R_4*P(:,4) - d2;
end