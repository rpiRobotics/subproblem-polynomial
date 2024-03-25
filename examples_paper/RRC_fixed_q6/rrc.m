kin = define_RRC_fixed_q6_rational();
kin.P(:,end) = 0;

%% Pose
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

R_06 = rot(ey, deg2rad(30))*rot(ez, deg2rad(20))*rot(kin.H(:,6), deg2rad(10));
R_06_num = double(R_06);
p_0T_num = [30.12;20.34;20.56];

%% Search-based IK
kin_num = kin;
kin_num.P = double(kin.P);
kin_num.H = double(kin.H);


[Q, is_LS_vec] = IK.IK_2_intersecting(R_06_num, p_0T_num, kin_num)
xlabel("q_4")
ylabel("e(q_4)")


[R_06_t, T_t] = fwdkin(kin, Q(:,1));

q = Q(:,1)
vpa(R_06_t - R_06)
vpa(T_t - p_0T)
tan(Q(4,:)/2)

%% Symbolic rational pose

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

syms xa xb xg
%R_06 = half_tan_rot(ey, xg)*half_tan_rot(ez, xb)*half_tan_rot(kin.H(:,6), xa);
R_06 = half_tan_rot(ey, xg)*half_tan_rot(ez, xb);

syms t1 t2 t3 real
p_0T = [t1 t2 t3]';
%% Rational approximations for tan half angle
[Nb, Db] = rat(tan(deg2rad(20)/2), 1e-6)
[Ng, Dg] = rat(tan(deg2rad(30)/2), 1e-6)

%%
polynomial_IK.two_intersecting(kin,R_06, p_0T, 'RRC_fixed_q6_eqns.txt')

%% Find remaining joint angles

x4_vec = [-1.1174935481027253756, -0.95143823714605743324, 0.96588325167689505419, 1.1915912313344271761, 1.7215238628117447496, 2.2154111526469432716];
q4_vec = 2*atan(x4_vec);

Q = [];
for i = 1:width(q4_vec)
    [e_vec, Q_partial] = alignment_err_given_q4(q4_vec(i), p_0T_num, R_06_num, kin_num);
    q_partial = Q_partial(:, abs(e_vec) == min(abs(e_vec)));
    Q(:,i) = q_from_partial(q_partial, kin_num, R_06_num);
end

Q

%% Double check forward kin
for i = 1:width(Q)
    [R_t, p_T] = fwdkin(kin_num, Q(:,1));
    [R_t p_T] - [R_06_num p_0T_num]
end

%%
function [e_vec, Q_partial] = alignment_err_given_q4(q4, p_16, R_06, kin)
    e_vec = NaN([1 4]);
    Q_partial = NaN([4 4]);
    % find up to 4 solutions of (q1, q2, q3) using Subproblem 5
    p_35_3 = kin.P(:,4) + rot(kin.H(:,4), q4)*kin.P(:,5);
    
    [t1, t2, t3] = subproblem.sp_5( ...
                   -kin.P(:,2), p_16, kin.P(:,3), p_35_3, ...
                   -kin.H(:,1), kin.H(:,2), kin.H(:,3));

    for i_q123 = 1:length(t1)
        q1 = t1(i_q123);
        q2 = t2(i_q123);
        q3 = t3(i_q123);

        R_04 = rot(kin.H(:,1),q1) * rot(kin.H(:,2),q2) ...
             * rot(kin.H(:,3),q3) * rot(kin.H(:,4),q4);

        e_i = kin.H(:,5)' * R_04' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);
        e_vec(i_q123) = e_i;
        Q_partial(:,i_q123) = [q1; q2; q3; q4];
    end
end

function q_i = q_from_partial(q_partial, kin, R_06)
    R_04 = rot(kin.H(:,1),q_partial(1)) * rot(kin.H(:,2),q_partial(2)) ...
         * rot(kin.H(:,3),q_partial(3)) * rot(kin.H(:,4),q_partial(4));
    [q5, q5_is_LS] = subproblem.sp_1(kin.H(:,6), R_04'*R_06*kin.H(:,6),  kin.H(:,5));
    [q6, q6_is_LS] = subproblem.sp_1(kin.H(:,5), R_06'*R_04*kin.H(:,5), -kin.H(:,6));
    q_i = [q_partial; q5; q6];
end
