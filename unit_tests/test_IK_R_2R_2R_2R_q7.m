kin = robot_kin.sawyer;
kin.P(:,end) = 0;
kin.P(:,1) = 0;

%%
R_07 = eye(3);
p_0T = [0.5; 0.5; 0.25];

SEW = sew_conv([0;0;1]);
psi = 0;
%%

q = [1 2 3 4 5 6 7]/10

[R_07, p_0T, T_SEW] = fwdkin_inter(kin, q, [1,4,7])
psi = SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3))


%%
tic
[Q, is_LS_vec] = IK_R_2R_2R_2R_q7(R_07, p_0T, SEW, psi, kin, true)
toc
%%
psi_vec = NaN(1,width(Q));
for i = 1:width(Q)
    [R_t, T_t, T_SEW] = fwdkin_inter(kin, Q(:,i), [1,4,7]);
    psi_vec(i) = SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3));

end
psi_vec

%%
[Q_t, is_LS_vec_t] = sawyer_find_remaining_angles(R_07, p_0T, SEW, psi, kin, Q(7,:))

