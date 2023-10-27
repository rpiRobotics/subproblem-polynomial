function [kin, q_given, T_given, R_given] = define_hyperplane_bot()
pi_sym = sym(pi);

alpha_vec = [-pi_sym/6 pi_sym/6 pi_sym/3 -pi_sym/6 -pi_sym/6 pi_sym/4];
a_vec     = [10 100 150 20 5 20];
d_vec     = [100 50 50 -50 -20 10];

q_given = [-pi_sym/6 pi_sym/2 -pi_sym/3 pi_sym/2 pi_sym/6 -pi_sym/6];

kin = dh_to_kin(alpha_vec, a_vec, d_vec);

% [R_06, T] = fwdkin(kin, q_given);
% T = simplify(T)
% R_06 = simplify(R_06)

T_given = [3395/128+7005*sqrt(3)/64
           4155*sqrt(3)/128 + 5455/64
           1845*sqrt(3)/64 - 695/32];

R_1 = [-59*sqrt(3)/512 + 177/256
       57*sqrt(3)/256 + 39/512
       15*sqrt(3)/128 + 137/256];
R_2 = [-sqrt(2)*(114*sqrt(3)+433)/1024
        sqrt(2)*(39* sqrt(3)+122)/1024
       -sqrt(2)*(55* sqrt(3)-246)/512];
R_3 = [ sqrt(2)*(66* sqrt(3)-115)/1024
        sqrt(2)*(21* sqrt(3)-650)/1024
        sqrt(2)*(59* sqrt(3)+90) /512];

% R = R_06 * kin.RT;
R_given = [R_1 R_2 R_3];
end