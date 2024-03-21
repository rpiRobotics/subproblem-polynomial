%% Radical definition

kin = define_RRC_fixed_q6();

kin.P(:,end) = 0;
kin.P(1,5) = round(kin.P(1,5)*1000)/1000;

%% Rational definition

kin = define_RRC_fixed_q6_rational();
kin.P(:,end) = 0;
%%
q = zeros([7 1]);

%% Symbolic pose
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

syms a b g
alpha = a;
beta = b;
gamma = g;
R_06 = rot(ey, gamma)*rot(ez, beta)*rot(kin.H(:,6), alpha);

assert(dot(ez,kin.H(:,6)) == 0)

syms t1 t2 t3 real
p_0T = [t1 t2 t3]' + kin.P(:,1) - R_06*kin.P(:,7);


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

%% Easy pose

R_06 = eye(3);
p_0T = [30;20;20];
%% Hard pose

R_06 = rot(ey, deg2rad(30))*rot(ez, deg2rad(20))*rot(kin.H(:,6), deg2rad(10));
R_06 = double(R_06);
p_0T = [30.12;20.34;20.56];
%%
[Nb, Db] = rat(tan(deg2rad(20)/2), 1e-6)
[Ng, Dg] = rat(tan(deg2rad(30)/2), 1e-6)
%%
R_06 = rot(ey, deg2rad(0))*rot(ez, deg2rad(20))*rot(kin.H(:,6), deg2rad(10));
R_06 = double(R_06);
p_0T = [30.12;20.34;20.56];
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
%% Search-based IK
kin_num = kin;
kin_num.P = double(kin.P);
kin_num.H = double(kin.H);


[Q, is_LS_vec] = IK.IK_2_intersecting(R_06, p_0T, kin_num)
xlabel("q_4")
ylabel("e(q_4)")


[R_06_t, T_t] = fwdkin(kin, Q(:,1));

q = Q(:,1)
vpa(R_06_t - R_06)
vpa(T_t - p_0T)
tan(Q(4,:)/2)

%%
polynomial_IK.two_intersecting(kin,R_06, p_0T, 'RRC_fixed_q6_eqns.txt')






%%
syms theta
x4 = tan(theta/2);
P = x4^10 + 7.9957332906938004840*x4^9 + 71.829437902169774489*x4^8 - 320.39582634852254622*x4^7 - 143.45540174223805211*x4^6 + 1019.3814394286641758*x4^5 + 175.33406704389630542*x4^4 - 1133.1028305538653591*x4^3 - 271.64831216219044916*x4^2 + 424.58471379419422059*x4 + 161.60322241701215186
hold on
fplot(-1e-2*P, [-pi pi])
ylim([-1.5 1])
xlim([-pi pi])
hold off
% yline(0)

%%
syms x4
P = x4^10 + 7.9957332906938004840*x4^9 + 71.829437902169774489*x4^8 - 320.39582634852254622*x4^7 - 143.45540174223805211*x4^6 + 1019.3814394286641758*x4^5 + 175.33406704389630542*x4^4 - 1133.1028305538653591*x4^3 - 271.64831216219044916*x4^2 + 424.58471379419422059*x4 + 161.60322241701215186
xsolns = roots(sym2poly(P))
complex_solns = xsolns(imag(xsolns) ~= 0)
xline(2*atan(real(complex_solns)), 'r', LineWidth=2)