[kin, q_given, T_given, R_given] = define_hyperplane_bot();

[R_06, p_0T] = fwdkin(kin, q_given);
% R_0T = R_06 * kin.RT;kin
R_06 = simplify(R_06)
p_0T = simplify(p_0T)

%% Try home position s.t. h_5 = ez
% [theta1, theta2, is_LS] = subproblem.sp_2([0;0;1], kin.H(:,5), -kin.H(:,4), kin.H(:,3))
kin.H(:,[3 4 ]) = kin.H(:,[4 3])

%% Try using rational numbers for everything

kin.P = round(kin.P*1000)/1000
kin.H = round(kin.H*1000)/1000
p_0T = round(p_0T*1000)/1000
R_06 = round(R_06*1000)/1000


%% Plot robot pose
diagrams.setup;
hold on
options = {"show_arrows", true,...
           "unit_size", 20,...
           "show_arrow_labels", false,...
           "show_joint_labels", false,...
           "show_base_label", false,...
           "show_base_frame", false,...
           "show_task_frame", false,...
           "cyl_half_length", 10,...
           "cyl_radius", 5};

diagrams.robot_plot(kin, q_given, options{:});

hold off
diagrams.redraw;

%% Make sure search-based IK works
kin_double.joint_type = kin.joint_type;
kin_double.H = double(kin.H);
kin_double.P = double(kin.P);

[Q, is_LS_vec] = IK.IK_gen_6_dof_mex(double(R_06), double(p_0T), kin_double)


[R_06_t, T_t] = fwdkin(kin_double, Q(:,1))
%%
R_06 = eye(3);
% T = [300; 70; 160];
p_0T = [20 200 10]';
%% Approx kinematics and pose

kin.P = round(kin.P, 6, "significant");
kin.H = round(kin.H, 6, "significant");
kin.RT = round(kin.RT, 6, "significant");

R_06 = round(R_06, 6, "significant");
p_0T = round(p_0T, 6, "significant");

%%
kin.P = double(kin.P);
kin.H = double(kin.H);
kin.RT = double(kin.RT);

R_06 = double(R_06);
p_0T = double(p_0T);

%%
polynomial_IK.general_6R(kin, R_06, p_0T, 'hyperplane_eqns.txt');

%% Alternate formulation: search and q1 qnd q6
p_0T = p_0T;
p_06 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);

syms x1 x2 x4 x6 real

R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_34 = half_tan_rot(kin.H(:,4), x4);
R_56 = half_tan_rot(kin.H(:,6), x6);

[lhs_1, lhs_2, R_23] = half_tan_sp5_R( ...
    -kin.P(:,3), ...
    R_01'*(p_06 - R_06*R_56'*kin.P(:,6))-kin.P(:,2), ...
    kin.P(:,4), ...
    kin.P(:,5), ...
    -kin.H(:,2), kin.H(:,3), kin.H(:,4), [x2 x4]);

% lhs_err_1 = kin.H(:,1)' * R_12 * R_23 * R_34 * kin.H(:,5) - kin.H(:,1)' * R_06 * R_56' * kin.H(:,5); % no x1
lhs_err_2 = kin.H(:,2)' * R_23 * R_34 * kin.H(:,5) - kin.H(:,2)' * R_01' * R_06 * R_56' * kin.H(:,5); % no x2
lhs_err_3 = kin.H(:,3)' * R_34 * kin.H(:,5) - kin.H(:,3)' * (R_01 * R_12)' * R_06 * R_56' * kin.H(:,5); % no x3
% lhs_err_4 = kin.H(:,4)' * kin.H(:,5) - kin.H(:,4)' * (R_01 * R_12 * R_23)' * R_06 * R_56' * kin.H(:,5); % no x4


eqns_lhs = [lhs_err_2; lhs_err_3; lhs_1; lhs_2];

eqns_lhs = simplify(eqns_lhs)

eqns_frac = simplifyFraction(eqns_lhs);
eqns_num = numden(eqns_frac);

fileID = fopen('hyperplane_eqns_alt1.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    fprintf(fileID, string(eqns_num(i)) + '\n');
end

%% Alternate formulation: search and q1 and q2, but sp5 on (q4,q5,q6)
p_0T = p_0T;
p_06 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);

syms x1 x2 x4 x6 real

R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_34 = half_tan_rot(kin.H(:,4), x4);
R_56 = half_tan_rot(kin.H(:,6), x6);

[lhs_1, lhs_2, R_23] = half_tan_sp5_R( ...
    -kin.P(:,3), ...
    R_01'*(p_06 - R_06*R_56'*kin.P(:,6))-kin.P(:,2), ...
    kin.P(:,4), ...
    kin.P(:,5), ...
    -kin.H(:,2), kin.H(:,3), kin.H(:,4), [x2 x4]);

% lhs_err_1 = kin.H(:,1)' * R_12 * R_23 * R_34 * kin.H(:,5) - kin.H(:,1)' * R_06 * R_56' * kin.H(:,5); % no x1
lhs_err_2 = kin.H(:,2)' * R_23 * R_34 * kin.H(:,5) - kin.H(:,2)' * R_01' * R_06 * R_56' * kin.H(:,5); % no x2
lhs_err_3 = kin.H(:,3)' * R_34 * kin.H(:,5) - kin.H(:,3)' * (R_01 * R_12)' * R_06 * R_56' * kin.H(:,5); % no x3
% lhs_err_4 = kin.H(:,4)' * kin.H(:,5) - kin.H(:,4)' * (R_01 * R_12 * R_23)' * R_06 * R_56' * kin.H(:,5); % no x4


eqns_lhs = [lhs_err_2; lhs_err_3; lhs_1; lhs_2];

eqns_lhs = simplify(eqns_lhs)

eqns_frac = simplifyFraction(eqns_lhs);
eqns_num = numden(eqns_frac);

fileID = fopen('hyperplane_eqns_alt2.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    fprintf(fileID, string(eqns_num(i)) + '\n');
end
