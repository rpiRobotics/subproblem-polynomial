kin = robot_kin.sawyer; % From stereo-sew repo

kin.P(:,end) = 0;
kin.P(:,1) = 0;

kin.P = sym(kin.P);
kin.H = sym(kin.H);

%% Define end effector pose
R_07 = eye(3);
p_0T = [0.5; 0.5; 0.25];

SEW = sew_conv([0;0;1]);
psi = 0;

%% Search-based (Searching over theta_w)
kin_num = robot_kin.sawyer;
kin_num.P(:,end) = 0;
kin_num.P(:,1) = 0;

[Q, is_LS_vec] = SEW_IK.IK_R_2R_2R_2R(R_07, p_0T, SEW, psi, kin_num, true)

% Expected values for q_6 and x_7
sort(Q(7,:))
sort(tan(Q(7,:)/2))

q_num = Q(:,1);
%% Search-based (Searching over q_7)
[Q, is_LS_vec] = IK_R_2R_2R_2R_q7(R_07, p_0T, SEW, psi, kin_num, true)

% Expected values for q_6 and x_7
% (including spurious solutions with elbow on wrong half plane)
sort(Q(7,:))
sort(tan(Q(7,:)/2))
%% Check a given solution to see if forward kinematics matches up

i_test = 6

[R_t, T_t, T_SEW] = fwdkin_inter(kin_num, Q(:,i_test), [1,4,7])
SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3))

%% Find the system of 4 multivariate polynomials

syms x1 x2 x6 x7 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_56 = half_tan_rot(kin.H(:,6), x6);
R_67 = half_tan_rot(kin.H(:,7), x7);

% q6 (SP4 based on q7 search)
p_07 = p_0T - R_07 * kin.P(:,8);
[~, n_SEW] = SEW.inv_kin([0;0;0], p_07, psi);

[~, a, b, c] = half_tan_sp4(R_67*R_07'*n_SEW, kin.P(:,6), -kin.H(:,6), n_SEW'*p_07);
lhs_x6 = a*x6^2 + b*x6 + c;

p_WE = R_07 * R_67' * R_56' * kin.P(:,6);

% q1
[~, a, b, c] = half_tan_sp3(kin.P(:,2), p_07-p_WE, kin.H(:,1), norm(kin.P(:,4)));
lhs_x1 = a*x1^2 + b*x1 + c;

% q2 (SP 4)
[~, a, b, c] = half_tan_sp4(kin.H(:,3), R_01'*(p_07-p_WE)-kin.P(:,2), -kin.H(:,2), kin.H(:,3)'*kin.P(:,4));
lhs_x2 = a*x2^2 + b*x2 + c;

% q3 (SP 1)
R_23 = half_tan_sp1_R(kin.P(:,4), R_12' * (R_01' * (p_07-p_WE) - kin.P(:,2)), kin.H(:,3));

% Error eqn
lhs_err = kin.H(:,4)' * (R_01 * R_12 * R_23)' * R_07 * (R_56 * R_67)'* kin.H(:,5) - kin.H(:,4)' * kin.H(:,5);


eqns_lhs = [lhs_x1; lhs_x2; lhs_x6; lhs_err];
eqns_lhs = simplify(eqns_lhs)

eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)
%% Export as a text file

filename = "sawyer_polynomials.txt";

fileID = fopen(filename,'w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i := ", i);
    % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + ':\n\n');
end

%% Sanity check: Make sure eqns are correct by checking with 1D search results
x_num = tan(q_num([1 2 6 7])'/2)
vpa(subs(eqns_lhs, [x1 x2 x6 x7], x_num))
vpa(subs(eqns_frac, [x1 x2 x6 x7], x_num))
vpa(subs(eqns_num, [x1 x2 x6 x7], x_num))


%% Given results from Maple, find remaining angles
%  (and filter out spurious solutions)

x7 = [-3.4051036718467858336, -2.3561794409947539619, -2.2437799661048489284, -1.3568363548342626941, -1.1844834747020878184, -0.42917381786918497380, -0.36743065645254205841, -0.29882092870136799522, -0.056161861197657397506, 0.043058763171974378376, 0.048163473698596596043, 0.29482858711463740550, 0.53113392608303378331, 0.84643436794217826784, 0.86263792375808204401, 1.1185307067448818295]
q7_vec = 2*atan(x7)

[Q, is_LS_vec] = sawyer_find_remaining_angles(R_07, p_0T, SEW, psi, kin_num, q7_vec)


psi_vec = NaN(1,width(Q));
for i = 1:width(Q)
    [R_t, T_t, T_SEW] = fwdkin_inter(kin_num, Q(:,i), [1,4,7]);
    psi_vec(i) = SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3));

end
psi_vec

Q_correct = Q(:,abs(psi_vec)<0.1) % Each q as a column
vpa(Q_correct', 10) % Each q as a row (like in paper)

%% Check all IK solutions
psi_vec = NaN(1,width(Q_correct));
for i = 1:width(Q_correct)
    [R_t, T_t, T_SEW] = fwdkin_inter(kin_num, Q_correct(:,i), [1,4,7])
    psi_vec(i) = SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3));
end
psi_vec
