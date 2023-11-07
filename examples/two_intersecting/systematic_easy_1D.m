kin = define_easy_1D_unit();
%%
R_06 = eye(3);
p_06 = [1;1;1]*1.25;

[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);

%%
kin = define_easy_1D_varying();

R_06 = eye(3);
p_06 = [1;1;1]*2;

[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);

%%
[R_t, P_t] = fwdkin(kin, Q(:,1))
%%
syms T1 T2 T3 real
syms R [3 3] real

R_06 = R
p_06 = [T1 T2 T3]'
%%
polynomial_IK.two_intersecting(kin, R_06, p_06, 'easy_1D_eqns.txt')

%%
[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_06, kin);
s = [-2.37742, -0.773904, -0.343369, 0.0691691, 0.198726, 1.74226]
xline(2*atan(s), 'r', LineWidth=2);