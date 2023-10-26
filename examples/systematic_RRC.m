zv = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

% Define in inches,
p01 = zv;
p12 = 20*ex -4*ey;
p23 = 4*ey;
p34 = sym(21.5)*ex + sym(3.375)*ey;
p45 = -sym(3.375)*ey;
p56 = sym(21.5)*ex+sym(3.325)*ey;
p67 = zv;

kin.P = [p01 p12 p23 p34 p45 p56 p67];
kin.H = [ex ez ex ez ex ez];
kin.joint_type=[0 0 0 0 0 0];
%%
q = zeros([7 1]);
%% Plot robot pose
diagrams.setup;
hold on
options = {"show_arrows", true,...
           "unit_size", 5,...
           "show_arrow_labels", false,...
           "show_joint_labels", false,...
           "show_base_label", false,...
           "show_base_frame", false,...
           "show_task_frame", false,...
           "cyl_half_length", 2,...
           "cyl_radius", 1};

diagrams.robot_plot(kin, q, options{:});

hold off
diagrams.redraw;
%%

kin_num = kin;
kin_num.P = double(kin.P);

R_06 = eye(3);
p_0T = [20;20;20];

[Q, is_LS_vec] = IK.IK_gen_6_dof(R_06, p_0T, kin_num)


[R_06_t, T_t] = fwdkin(kin, Q(:,1))

q = Q(:,1)

%%
p_06 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);


%%%%%%%%%%%%
% p_06 = round(p_06, 6, "significant");
% % R_06 = round(R_06, 6, "significant");
% kin.H = round(kin.H, 6, "significant");
% kin.P = round(kin.P, 6, "significant");
%%%%%%%%%%%%

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
%%
eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)
%%
fileID = fopen('RRC_eqns.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    % fprintf(fileID, string(vpa(eqns_num(i), 6)) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end
