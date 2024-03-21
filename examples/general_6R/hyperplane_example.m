%%
% Tolerances for first working attempt:
% beta / gamma : 1e-3
% p_0T: 1e-4
% H: 1e-6
% P: 1e-6

% Tolerances for more accurate attempt:
% All 1e-6
%%
[kin, q_given, T_given, R_given] = define_hyperplane_bot();

[R_06, p_0T] = fwdkin(kin, q_given);
% R_0T = R_06 * kin.RT;kin
R_06 = simplify(R_06)
p_0T = simplify(p_0T)

%% Cut out last part of P
p_0T = p_0T - R_given*kin.RT'*kin.P(:,end);
kin.P(:,end) = 0;

%% Parameterize R_06 and approximate p_0T
% R_06 = R_z R_y R_z
%      = R_gamma R_beta R_alpha
% R_z^T R_06 e_z = R_y ez 
R_06 = double(R_given*kin.RT');

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
[gamma, beta] =  subproblem.sp_2(R_06 * ez, ez, -ez, ey)

R_06 * ez - rot(ez,gamma(1)) * rot(ey, beta(1)) * ez

x_beta_exact  = tan(beta(1)/2);
x_gamma_exact = tan(gamma(1)/2);
[N_xb, D_xb] = rat(x_beta_exact, 1e-6);
[N_xg, D_xg] = rat(x_gamma_exact, 1e-6);
x_beta = sym(N_xb) / sym(D_xb)
x_gamma = sym(N_xg) / sym(D_xg)

R_06_rat = half_tan_rot(ez, x_gamma) * half_tan_rot(ey, x_beta)
double(R_06_rat * ez - R_06 * ez)

% Also find R_alpha precisely
alpha = subproblem.sp_1(ey, R_06_rat'*R_06*ey, ez)
R_alpha = rot(ez, alpha)
double((R_06_rat * R_alpha)'*R_06-eye(3))

double(p_0T)
[N_P, D_P] = rat(p_0T, 1e-4);
p_0T = N_P./D_P
double(p_0T)
%% Rational Approximation for kinematics
% exact H is as follows:
% s = sin(pi/6) = 1/2
% c = cos(pi/6) = sqrt(3)/2
% 
% H = 
% [0  0  0  0  0  0
%  0  s  0 -c -s  0
%  1  c  1  s  c  1]
% Use an angle slightly off from pi/6 instead

x_exact = tan((sym(pi)/6)/2)
[N, D] = rat(x_exact, 1e-6) % TODO
x = N/D
s = 2*x/(x^2+1)
c = (-x^2+1)/(x^2+1)

e_s = double(s - sin(pi/6))
e_c = double(c - cos(pi/6))

s^2+c^2

kin.H = ...
[0  0  0  0  0  0
 0  s  0 -c -s  0
 1  c  1  s  c  1];

[NP, ND] = rat(kin.P, 1e-6); % TODO
kin.P = NP./ND;


%% Test out rational kinematics approx accuracy
[R_T, T_t] = fwdkin(kin, q_given);


double(R_T' * R_given*kin.RT' - eye(3))
double(R_T' * R_06_rat * R_alpha - eye(3))
double(T_t - p_0T )
% double(T_t - p_0T +R_given*kin.RT'*[20;0;10])

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

% [Q, is_LS_vec] = IK.IK_gen_6_dof_mex(double(R_06), double(p_0T), kin_double)
[Q, is_LS_vec] = IK.IK_gen_6_dof_mex(double(R_06_rat * R_alpha), double(p_0T), kin_double)

[R_06_t, T_t] = fwdkin(kin_double, Q(:,1))
%%
q1_list = uniquetol(sort(Q(1,:)), 1e-4)
[tan(q1_list/2 - 1e-3) 
    tan(q1_list/2)
    tan(q1_list/2 + 1e-3)]

%%
q2_list = uniquetol(sort(Q(2,:)), 1e-4)
[tan(q2_list/2 - 1e-3)
    tan(q2_list/2)
    tan(q2_list/2 + 1e-3)]
%%
% close
R_06 = eye(3);
% T = [300; 70; 160];
% p_0T = [30 150 10]';
% p_0T = round(rand_vec*200)
 p_0T=[   71-20
   103
    97-10];
%%
kin.P = double(kin.P);
kin.H = double(kin.H);
kin.RT = double(kin.RT);

R_06 = double(R_06);
p_0T = double(p_0T);

%%
% [eqns, eqns2] = polynomial_IK.general_6R(kin, R_06, p_0T, 'hyperplane_eqns.txt');
eqns_num = polynomial_IK.general_6R(kin, R_06_rat, p_0T, 'hyperplane_eqns.txt');

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
