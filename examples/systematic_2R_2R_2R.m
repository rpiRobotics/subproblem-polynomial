%% First a general 2R-2R-2R robot
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];
kin.H = rand_normal_vec(6);
kin.P = [rand_vec zv rand_vec zv rand_vec zv rand_vec];
kin.joint_type = zeros(1,6);

q = rand_angle([6 1]);

[R_06,p_0T] = fwdkin(kin, q)



%% zero offsets 

kin.H = [ez ey ez ey ez ey];
kin.P = [ez zv ex zv ex zv ex];
kin.joint_type = zeros(1,6);

R_06 = eye(3);
p_0T = [0; -1.5; 1.5];
%%
[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)


%%

% x3 in terms of x4
syms x4 real


p_16 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);
R_34 = half_tan_rot(kin.H(:,4), x4);
x3 = half_tan_sp3(R_34*kin.P(:,5), -kin.P(:,3), kin.H(:,3), norm(p_16))

% x1,x2 using subproblem 2
R_23 = half_tan_rot(kin.H(:,3), x3);
[x1, x2] = half_tan_sp2(p_16, kin.P(:,3)+R_23*R_34*kin.P(:,5), -kin.H(:,1), kin.H(:,2));

R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_04 = R_01 * R_12 * R_23 * R_34;
e = kin.H(:,5)' * R_04.' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);