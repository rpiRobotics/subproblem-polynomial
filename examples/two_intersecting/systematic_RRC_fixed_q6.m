kin = define_RRC_fixed_q6()
kin.P(:,end) = 0;
kin.P(1,5) = round(kin.P(1,5)*1000)/1000;
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
kin_num.H = double(kin.H);

R_06 = eye(3);
p_0T = [30;20;20];

[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_0T, kin_num)
xlabel("q_4")
ylabel("e(q_4)")


[R_06_t, T_t] = fwdkin(kin, Q(:,1));

q = Q(:,1)
vpa(R_06_t - R_06)
vpa(T_t - p_0T)
%%
polynomial_IK.two_intersecting(kin,R_06, p_0T, 'RRC_fixed_q6_eqns.txt')