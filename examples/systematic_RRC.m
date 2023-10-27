kin = define_RRC();
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
polynomial_IK.general_6R(kin, R_06, p_0T, 'RRC_eqns.txt');