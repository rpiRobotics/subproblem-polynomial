kin = robot_kin.sawyer

kin.P(:,end) = 0;
%% Easy pose
R_07 = eye(3);
p_0T = [0.25; 0.25; 0.25];

SEW = sew_conv([0;0;1]);
psi = 0;

%% Search-based

[Q, is_LS_vec] = SEW_IK.IK_R_2R_2R_2R(R_07, p_0T, sew, psi, kin, true)

%%

RAT_TOL = 1e-3;

syms x1 x2 x3 x4 x5 x_wrist real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_34 = half_tan_rot(kin.H(:,4), x4);
R_45 = half_tan_rot(kin.H(:,5), x5);

R_05 = R_01 * R_12 * R_23 * R_34 * R_45;

% Error eqn
lhs_err = kin.H(:,6)' * R_05' * R_07 * kin.H(:,7) - kin.H(:,6)' * kin.H(:,7);

% q1
p_W_EE_0 = kin.P(:,8);
W = p_0T - R_07 * p_W_EE_0;
S = kin.P(:,1);
p_17 = W-S;
h_1 = kin.H(:,1);

w_hat = vec_normalize(W - S);
[N, D] = rat(w_hat, RAT_TOL)
w_hat_rat = sym(N)./sym(D);

p_WE_0 = -kin.P(:,6);
d_WE = norm(p_WE_0);

p_E2_0 = -kin.P(:,4);
d_2E = norm(p_E2_0);

[~, n_SEW] = SEW.inv_kin(S, W, psi);
[N, D] = rat(n_SEW, RAT_TOL)
n_SEW_rat = sym(N)./sym(D);

[N, D] = rat(d_WE, RAT_TOL)
d_WE_rat = sym(N)./sym(D);

p_WE = half_tan_rot(n_SEW_rat, x_wrist) * (-w_hat_rat) * d_WE_rat; % TODO - make rational in a more clever way?

[~, a, b, c] = half_tan_sp3(kin.P(:,2), p_17+p_WE, h_1, d_2E);
lhs_x1 = a*x1^2 + b*x1 + c;

% (q2, q3)
R_10 = R_01';
R_21 = R_12';
R_32 = R_23';
[~, ~, a3, b3, c3, a2, b2, c2] = half_tan_sp2(kin.P(:,4), R_10*p_17+R_10*p_WE-kin.P(:,2), kin.H(:,3), -kin.H(:,2))
lhs_x2 = a2*x2^2 + b2*x2 + c2;
lhs_x3 = a3*x3^2 + b3*x3 + c3;

% (q4, q5)
[~, ~, a5, b5, c5, a4, b4, c4] = half_tan_sp2(kin.P(:,6), R_32*R_21*(R_10*p_17-kin.P(:,2))-kin.P(:,4), kin.H(:,5), -kin.H(:,4));
lhs_x5 = a5*x5^2 + b5*x5 + c5;
lhs_x4 = a4*x4^2 + b4*x4 + c4;

%%
eqns_lhs = [lhs_x1; lhs_x2; lhs_x3; lhs_x4; lhs_x5; lhs_err];
eqns_lhs = simplify(eqns_lhs)

eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

filename = "sawyer_polynomials.txt";

fileID = fopen(filename,'w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end

function n = vec_normalize(vec)
    n =  vec / norm(vec);
end