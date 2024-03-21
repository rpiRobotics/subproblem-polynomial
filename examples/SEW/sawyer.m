kin = robot_kin.sawyer

kin.P(:,end) = 0;
kin.P(:,1) = 0;

kin.P = sym(kin.P)
kin.H = sym(kin.H)
%% Easy pose
R_07 = eye(3);
p_0T = [0.5; 0.5; 0.25];

SEW = sew_conv([0;0;1]);
psi = 0;

%% Search-based
% (But searching theta_w rather than q_7)
kin = robot_kin.sawyer

kin.P(:,end) = 0;
kin.P(:,1) = 0;
[Q, is_LS_vec] = SEW_IK.IK_R_2R_2R_2R(R_07, p_0T, SEW, psi, kin, true)
q_num = Q(:,2)

sort(tan(Q(7,:)/2))
%%
[R_t, T_t, T_SEW] = fwdkin_inter(kin, Q(:,5), [1,4,7])
SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3))

%% OLD APPROACH
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
%%
eqns_lhs = [lhs_x1; lhs_x2; lhs_x6; lhs_err];
eqns_lhs = simplify(eqns_lhs)

eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

filename = "sawyer_polynomials.txt";

fileID = fopen(filename,'w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i := ", i);
    % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + ':\n\n');
end

%% NEW APROACH
% syms x1 x2 x3 x6 x7 real
% R_01 = half_tan_rot(kin.H(:,1), x1);
% R_12 = half_tan_rot(kin.H(:,2), x2);
% R_23 = half_tan_rot(kin.H(:,3), x3);
% R_56 = half_tan_rot(kin.H(:,6), x6);
% R_67 = half_tan_rot(kin.H(:,7), x7);
% 
% p_07 = p_0T - R_07 * kin.P(:,8);
% [~, n_SEW] = SEW.inv_kin([0;0;0], p_07, psi);
% 
% % q6 (SP4 based on q7 search)
% [~, a, b, c] = half_tan_sp4(R_67*R_07'*n_SEW, kin.P(:,6), -kin.H(:,6), n_SEW'*p_07);
% lhs_x6 = a*x6^2 + b*x6 + c;
% 
% p_WE = R_07 * R_67' * R_56' * kin.P(:,6);
% 
% % q1
% [~, a, b, c] = half_tan_sp3(kin.P(:,2), p_07-p_WE, kin.H(:,1), norm(kin.P(:,4)));
% lhs_x1 = a*x1^2 + b*x1 + c;
% 
% % Find (q_2, q_3) with Subproblem 2
% [~, ~, a3, b3, c3, a2, b2, c2] = half_tan_sp2(kin.P(:,4), R_01'*p_07-R_01'*p_WE-kin.P(:,2), kin.H(:,3), -kin.H(:,2));
% lhs_x3 = a3*x3^2 + b3*x3 + c3;
% lhs_x2 = a2*x2^2 + b2*x2 + c2;
% 
% % Error eqn
% lhs_err = kin.H(:,4)' * (R_01 * R_12 * R_23)' * R_07 * (R_56 * R_67)'* kin.H(:,5) - kin.H(:,4)' * kin.H(:,5);

%%
% eqns_lhs = [lhs_x1; lhs_x2; lhs_x3; lhs_x6; lhs_err];
% eqns_lhs = simplify(eqns_lhs)
% 
% eqns_frac = simplifyFraction(eqns_lhs)
% eqns_num = numden(eqns_frac)
% 
% filename = "sawyer_polynomials_alt.txt";
% 
% fileID = fopen(filename,'w');
% for i = 1:length(eqns_num)
%     fprintf(fileID, "p%i := ", i);
%     % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
%     fprintf(fileID, string(eqns_num(i)) + ':\n\n');
% end




%% Make sure eqns are correct
x_num = tan(q_num([1 2 6 7])'/2)
vpa(subs(eqns_lhs, [x1 x2 x6 x7], x_num))
vpa(subs(eqns_frac, [x1 x2 x6 x7], x_num))
vpa(subs(eqns_num, [x1 x2 x6 x7], x_num))

%%
R_67 = rot(kin.H(:,7), q_num(7));

% verify q6
t6 = subproblem.sp_4(R_67*R_07*n_SEW_rat, kin.P(:,6), -kin.H(:,6), n_SEW_rat'*kin.P(:,7));
t6 = double(t4)
 q_num(6)
%%
% verify q1
                    %(kin.P(:,2),       p_17-p_WE,                      kin.H(:,1),    norm_p34)
t1 = subproblem.sp_3(kin.P(:,2), subs(p_17-p_WE, [x1 x2 x6 x7], x_num), kin.H(:,1), norm(kin.P(:,4)));
vpa(t1)
vpa(q_num(1))
%%
% verify q2
t2 = subproblem.sp_4(kin.H(:,3), subs(R_01'*(p_17-p_WE)-kin.P(:,2), [x1 x2 x6 x7], x_num), -kin.H(:,2), kin.H(:,3)'*kin.P(:,4));
double(t2)
q_num(2)
%%
% verify q3
t3 = subproblem.sp_1(kin.P(:,4), R_12' * (R_01' * (p_17-p_WE) - kin.P(:,2)), kin.H(:,3))
double(subs(t3, [x1 x2 x6 x7], x_num))
q_num(3)

% verify err
vpa(subs(lhs_err, [x1 x2 x6 x7], x_num))

%% Given results from Maple, find remaining angles

x7 = [-3.4051036718467858336, -2.3561794409947539619, -2.2437799661048489284, -1.3568363548342626941, -1.1844834747020878184, -0.42917381786918497380, -0.36743065645254205841, -0.29882092870136799522, -0.056161861197657397506, 0.043058763171974378376, 0.048163473698596596043, 0.29482858711463740550, 0.53113392608303378331, 0.84643436794217826784, 0.86263792375808204401, 1.1185307067448818295]
q7_vec = 2*atan(x7)

[Q, is_LS_vec] = sawyer_find_remaining_angles(R_07, p_0T, SEW, psi, kin, q7_vec)


psi_vec = NaN(1,width(Q));
for i = 1:width(Q)
    [R_t, T_t, T_SEW] = fwdkin_inter(kin, Q(:,i), [1,4,7]);
    psi_vec(i) = SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3));

end
psi_vec

Q_correct = Q(:,abs(psi_vec)<0.1)
vpa(Q_correct', 10)

%%
for i = 1:width(Q_correct)
    [R_t, T_t, T_SEW] = fwdkin_inter(kin, Q_correct(:,i), [1,4,7])
    psi_vec(i) = SEW.fwd_kin(T_SEW(:,1),T_SEW(:,2),T_SEW(:,3));
end
psi_vec

%%
function n = vec_normalize(vec)
    n =  vec / norm(vec);
end