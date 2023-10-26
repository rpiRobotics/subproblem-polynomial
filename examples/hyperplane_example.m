pi_sym = sym(pi);

alpha_vec = [-pi_sym/6 pi_sym/6 pi_sym/3 -pi_sym/6 -pi_sym/6 pi_sym/4];
a_vec     = [10 100 150 20 5 20];
d_vec     = [100 50 50 -50 -20 10];

q_given = [-pi_sym/6 pi_sym/2 -pi_sym/3 pi_sym/2 pi_sym/6 -pi_sym/6];

kin = dh_to_kin(alpha_vec, a_vec, d_vec);

[R_06, T] = fwdkin(kin, q_given);
T = simplify(T)
R_06 = simplify(R_06)

T_given = [3395/128+7005*sqrt(3)/64
           4155*sqrt(3)/128 + 5455/64
           1845*sqrt(3)/64 - 695/32]

R_1 = [-59*sqrt(3)/512 + 177/256
       57*sqrt(3)/256 + 39/512
       15*sqrt(3)/128 + 137/256];
R_2 = [-sqrt(2)*(114*sqrt(3)+433)/1024
        sqrt(2)*(39* sqrt(3)+122)/1024
       -sqrt(2)*(55* sqrt(3)-246)/512];
R_3 = [ sqrt(2)*(66* sqrt(3)-115)/1024
        sqrt(2)*(21* sqrt(3)-650)/1024
        sqrt(2)*(59* sqrt(3)+90) /512];

R = R_06 * kin.RT
R_given = [R_1 R_2 R_3]

%% Try using rational numbers for everything

kin.P = round(kin.P*1000)/1000
kin.H = round(kin.H*1000)/1000
T = round(T*1000)/1000
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

[Q, is_LS_vec] = IK.IK_gen_6_dof_mex(double(R_06), double(T), kin_double)


[R_06_t, T_t] = fwdkin(kin_double, Q(:,1))
%%
R_06 = eye(3);
% T = [300; 70; 160];
T = [20 200 10]';
%% Approx kinematics and pose

kin.P = round(kin.P, 6, "significant");
kin.H = round(kin.H, 6, "significant");
kin.RT = round(kin.RT, 6, "significant");

R_06 = round(R_06, 6, "significant");
T = round(T, 6, "significant");

%%
kin.P = double(kin.P);
kin.H = double(kin.H);
kin.RT = double(kin.RT);

R_06 = double(R_06);
T = double(T);
%%
p_0T = T;
p_06 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);

syms x1 x2 x3 x5 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
% R_34 = half_tan_rot(kin.H(:,4), x4);
R_45 = half_tan_rot(kin.H(:,5), x5);

[lhs_1, lhs_2, R_34] = half_tan_sp5_R( ...
    -kin.P(:,4), ...
    R_12'*(R_01'*p_06-kin.P(:,2))-kin.P(:,3), ...
    kin.P(:,5), ...
    kin.P(:,6), ...
    -kin.H(:,3), kin.H(:,4), kin.H(:,5), [x3 x5]);

% Error eqns
% R_05 = R_01 * R_12 * R_23 * R_34 * R_45;
%lhs_err = R_05 * kin.H(:,6) - R_06 * kin.H(:,6); % 3 equations
%lhs_err = (R_06 * kin.H(:,6))'*  R_05 * kin.H(:,6) -1
lhs_err_1 = kin.H(:,3)'*R_34*R_45*kin.H(:,6) - kin.H(:,3)'*(R_01*R_12)'*R_06*kin.H(:,6); % No x3
lhs_err_2 = kin.H(:,4)'*R_45*kin.H(:,6) - kin.H(:,4)'*(R_01*R_12*R_23)'*R_06*kin.H(:,6); % No x4

% lhs_err_3 = kin.H(:,1)'*R_12 * R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,1)'*R_06*kin.H(:,6); % no x1
% lhs_err_4 = kin.H(:,2)'*R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,2)'* R_01' * R_06*kin.H(:,6); % no x2
% lhs_err_5 = kin.H(:,5)'*kin.H(:,6) - kin.H(:,5)'*(R_01*R_12*R_23*R_34)'*R_06*kin.H(:,6); % No x5


eqns_lhs = [lhs_err_1; lhs_err_2; lhs_1; lhs_2];

eqns_lhs = simplify(eqns_lhs)
%% Alternate formulation: search and q1 qnd q6
p_0T = T;
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
%%
eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)
%%
fileID = fopen('hyperplane_eqns.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    % fprintf(fileID, string(vpa(eqns_num(i), 6)) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end

%% Test eqns

q = double(q_given);
x = tan(q/2)

double(subs(eqns_lhs(2), [x1 x2 x4 x6], x([1 2 4 6])))
double(subs(eqns_num(2), [x1 x2 x4 x6], x([1 2 4 6])))

%% Test eqns

q = double(q_given);
x = tan(q/2)

double(subs(eqns_lhs(4), [x1 x2 x3 x5], x([1 2 3 5])))
double(subs(eqns_num(4), [x1 x2 x3 x5], x([1 2 3 5])))

%%
R_01 = rot(kin.H(:,1), q_given(1));
R_12 = rot(kin.H(:,2), q_given(2));
R_23 = rot(kin.H(:,3), q_given(3));
R_34 = rot(kin.H(:,4), q_given(4));
R_56 = rot(kin.H(:,6), q_given(6));

l = double(-kin.P(:,3) + R_12'*(R_01'*(p_06 - R_06*R_56'*kin.P(:,6))-kin.P(:,2)))
r = double(R_23*kin.P(:,4) + R_23*R_34*kin.P(:,5))