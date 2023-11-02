ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin.joint_type = [0 0 0 0 0 0];
kin.P = [zv ex ex ex ex ex zv];
kin.H = [ez ey ez ey ez ey];

R_06 = eye(3);
p_06 = [1;1;1]*1.25;

[Q, is_LS_vec] = IK.IK_gen_6_dof(R_06, p_06, kin);
%%
[R_t, p_t] = fwdkin(kin, Q(:,2))
%%

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
R_05 = R_01 * R_12 * R_23 * R_34 * R_45;
%lhs_err = R_05 * kin.H(:,6) - R_06 * kin.H(:,6); % 3 equations
%lhs_err = (R_06 * kin.H(:,6))'*  R_05 * kin.H(:,6) -1
lhs_err_1 = kin.H(:,3)'*R_34*R_45*kin.H(:,6) - kin.H(:,3)'*(R_01*R_12)'*R_06*kin.H(:,6); % No x3
lhs_err_2 = kin.H(:,4)'*R_45*kin.H(:,6) - kin.H(:,4)'*(R_01*R_12*R_23)'*R_06*kin.H(:,6); % No x4

lhs_err_3 = kin.H(:,1)'*R_12 * R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,1)'*R_06*kin.H(:,6); % no x1
lhs_err_4 = kin.H(:,2)'*R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,2)'* R_01' * R_06*kin.H(:,6); % no x2
lhs_err_5 = kin.H(:,5)'*kin.H(:,6) - kin.H(:,5)'*(R_01*R_12*R_23*R_34)'*R_06*kin.H(:,6); % No x5


eqns_lhs = [lhs_err_1; lhs_err_2; lhs_err_3; lhs_err_4; lhs_err_5; lhs_1; lhs_2];
eqns_lhs = simplify(eqns_lhs)
%%
eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)

%% Test eqns

q = Q(:,1);
x = tan(q/2)

double(subs(eqns_lhs(5), [x1 x2 x3 x5], x([1 2 3 5])'))
double(subs(eqns_num(5), [x1 x2 x3 x5], x([1 2 3 5])'))

%%
fileID = fopen('easy_gen_6_dof_eqns.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    fprintf(fileID, string(vpa(eqns_num(i))) + '\n');
    % fprintf(fileID, string(eqns_num(i)) + '\n');
end

%%
x1s = [-1.2064,-1.1228,-0.859648,-0.431921,-0.289665,-0.221697,-0.119455,-0.0759565,0.0765156,0.078783,0.109625,0.110009,0.113658,0.113699,0.118011,0.11978,0.216534,0.299297,0.308924,0.335465,0.341205,0.377572,0.686284,0.78967,0.820456,0.820773,0.831854,0.834321,0.967919,0.975219,0.995004,1.00568,1.03023,1.17454,1.37358,1.4093,1.60625,1.63508,1.68159,1.83729,1.87804,2.10491,2.1127,2.11417,2.15083,2.18186,2.1889,2.19358,2.21036,2.25812,2.42502,2.52261,3.39834,3.44777]
q1s = 2*atan(x1s)