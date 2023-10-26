kin = define_CRX;
%R_06 = eye(3);
%p_0T = [0.25; 0.25; 0.25];
% [R_06, p_0T] = fwdkin(kin, rand_angle([7 1]))

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];
R_06 = rot(ez, deg2rad(30))*rot(ey, deg2rad(20))*rot(ex, deg2rad(10));
% R_06 = rot(ez, deg2rad(sym(30)))*rot(ey, deg2rad(sym(20)))*rot(ex, deg2rad(sym(10)));
% p_0T = [1415; 9265; 3589]/10000;
p_0T = [3141; 5926; 5358]/10000;


%%
% syms alpha beta gamma 
syms t1 t2 t3 real
%R_06 = rot(ez, gamma)*rot(ey, beta)*rot(ex, alpha);
syms R11 R12 R13 R21 R22 R23 R31 R32 R33  real
R_06 = [R11 R12 R13; R21 R22 R23; R31 R32 R33]
p_0T = [t1 t2 t3]' + kin.P(:,1) - R_06*kin.P(:,7);
%%
[Q, is_LS_vec] = IK_3_pairs_intersecting(double(R_06),p_0T,kin, true)

%%
[R_test, T_test] = fwdkin(kin, Q(:,1))
%%
syms x1 x2 x3 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_23 = half_tan_rot(kin.H(:,3), x3);
R_34 = half_tan_rot(kin.H(:,4), x4);

% Error eqn
R_04 = R_01 * R_12 * R_23 * R_34;
lhs_err = kin.H(:,5)' * R_04.' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);

% x3 eqn
p_16 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);
R_34 = half_tan_rot(kin.H(:,4), x4);
[~, a, b, c] = half_tan_sp3(R_34*kin.P(:,5), -kin.P(:,3), kin.H(:,3), norm(p_16));
lhs_x3 = a*x3^2 + b*x3 + c;

% x1 and x2 eqn
[~, ~, a1, b1, c1, a2, b2, c2] = half_tan_sp2(p_16, kin.P(:,3)+R_23*R_34*kin.P(:,5), -kin.H(:,1), kin.H(:,2));

lhs_x1 = a1*x1^2 + b1*x1 + c1;
lhs_x2 = a2*x2^2 + b2*x2 + c2;

eqns_lhs = [lhs_err; lhs_x1; lhs_x2; lhs_x3];
eqns_lhs = simplify(eqns_lhs)

%% Alternate approach with Subproblem 5
syms x1 x2 x4 real
R_01 = half_tan_rot(kin.H(:,1), x1);
R_12 = half_tan_rot(kin.H(:,2), x2);
R_34 = half_tan_rot(kin.H(:,4), x4);

p_16 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);
[lhs_1, lhs_2, R_23] = half_tan_sp5_R( ...
    -kin.P(:,3), ...
    R_01'*p_16, ...
    zv, ...
    kin.P(:,5), ...
    -kin.H(:,2), kin.H(:,3), kin.H(:,4), [x2 x4]);
R_04 = R_01 * R_12 * R_23 * R_34;
lhs_err = kin.H(:,5)' * R_04.' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);

eqns_lhs = [lhs_err; lhs_1; lhs_2];
eqns_lhs = simplify(eqns_lhs)

%%
eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

%%
fileID = fopen('crx_eqns.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    % fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end
%% Test eqns

q = Q(:,1);
x = tan(q/2);

double(subs(eqns_lhs(1), [x1 x2 x4], x([1 2 4])' ))
double(subs(eqns_num(1), [x1 x2 x4], x([1 2 4])' ))
%% Test eqns

q = Q(:,1);
x = tan(q/2);

double(subs(eqns_lhs(4), [x1 x2 x3 x4], x([1 2 3 4])' ))
double(subs(eqns_num(4), [x1 x2 x3 x4], x([1 2 3 4])' ))
%%
p1 = eqns_num(1)
p2 = eqns_num(2)
p3 = eqns_num(3)
p4 = eqns_num(4)

% dr = ResourceFunction["DixonResultant"][{p1, p2, p3, p4}, {x1, x2, x3}]
% Factor[dr]
% p_star = FactorList[dr][[12]][[1]]
% x4 /. NSolve[pstar[[1]], x4, Reals]
%%
x_roots = [-1.22679,-1.2064,-1.11032,-0.859648,-0.289665,-0.221697,-0.0587106,-0.0469804,0.0696124,0.0728794,0.109625,0.110009,0.113632,0.113664,0.118011,0.11978,0.216534,0.299297,0.308924,0.335465,0.341205,0.377572,0.78967,0.810497,0.820178,0.820773,0.831854,0.836245,0.888068,0.975219,0.995004,1.00933,1.03066,1.08293,1.17454,1.4093,1.68159,1.83729,1.87804,2.10455,2.1127,2.11546,2.15015,2.1889,2.19365,2.20063,2.21036,2.25812,2.37865,2.52261,3.41518,3.44777]
theta_roots = 2*atan(x_roots);

real_theta = theta_roots(imag(theta_roots)==0)

%%
[Q, is_LS_vec] = IK_3_pairs_intersecting(R_06,p_0T,kin, true)
% xline(real(theta_roots), 'r')
xline(real_theta, 'r')

%%
R_01 = rot(kin.H(:,1), q(1));
R_12 = rot(kin.H(:,2), q(2));
R_23 = rot(kin.H(:,3), q(3));
R_34 = rot(kin.H(:,4), q(4));

-kin.P(:,3) + R_12'*R_01'*p_16
R_23*R_34*kin.P(:,5)