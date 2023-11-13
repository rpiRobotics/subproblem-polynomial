%% Full 7-dof kinematics

kin_full = define_yumi();

%% 6-dof by fixing q_3 = pi/2
pi_sym = sym(pi);
[kin, R_6T] = fwdkin_partial(kin_full, pi_sym/2, 3);
%% Easy pose
R_06 = eye(3);

T =[50
 -355
  354];
%% Hard pose
R_06 = rot(ey, deg2rad(30))*rot(ez, deg2rad(20))*rot(kin.H(:,6), deg2rad(10));
R_06 = double(R_06);

T =[50.12
 -355.34
  354.56];
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
T = [t1 t2 t3]';

%% Symbolic rational pose

ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

syms xa xb xg
%R_06 = half_tan_rot(ey, xg)*half_tan_rot(ez, xb)*half_tan_rot(kin.H(:,6), xa);
R_06 = half_tan_rot(ey, xg)*half_tan_rot(ez, xb);

syms t1 t2 t3 real
T = [t1 t2 t3]';
%% Instead, fix q1
[kin, R_6T] = fwdkin_partial(kin_full, 0, 1);
kin.P(:,1) = 0;

%% Different q_3
[kin, R_6T] = fwdkin_partial(kin_full, pi_sym/6, 3);
kin.P = rot(ez,-pi_sym/6)*kin.P;
kin.H = rot(ez,-pi_sym/6)*kin.H;
%% Plot robot pose
diagrams.setup;
hold on
options = {"show_arrows", true,...
           "unit_size", 40,...
           "show_arrow_labels", false,...
           "show_joint_labels", false,...
           "show_base_label", false,...
           "show_base_frame", false,...
           "show_task_frame", false,...
           "cyl_half_length", 20,...
           "cyl_radius", 10};

% diagrams.robot_plot(kin_full, deg2rad([0 0 0 -90 180 0 0]), options{:});
% diagrams.robot_plot(kin_full, zeros([7 1]), options{:});
% diagrams.robot_plot(kin, deg2rad([0 0 -90 180 0 0]), options{:});
diagrams.robot_plot(kin_double, zeros([6 1]), options{:});

hold off
diagrams.redraw;
%% Make sure search-based IK works
kin_double.joint_type = kin.joint_type;
kin_double.H = double(kin.H);
kin_double.P = double(kin.P);

[Q, is_LS_vec] = IK.IK_gen_6_dof_mex(double(R_06), double(T), kin_double)


[R_06_t, T_t] = fwdkin(kin_double, Q(:,1))
view(0,-90)
pbaspect([1 1 1/5])
%%
q1_list = uniquetol(Q(1,:), 1e-2)
tan(q1_list/2)
%% Output system of multivariate polynomials
eqns = polynomial_IK.general_6R(kin, R_06, T, 'yumi_eqns.txt')

%% Test eqns
q = Q(:,1);
x = tan(q/2)

double(subs(eqns(4), [x1 x2 x3 x5], x([1 2 3 5])'))
double(subs(eqns(4), [x1 x2 x3 x5], x([1 2 3 5])'))

%% Get remaining joint angles
x_12 = [-3.183561363129973		0.43234343786735896
-2.505162523687448		0.3626464329484076
-0.4949608488236235		0.4617129265620604
-0.4188379662535908		0.607508293446981
 0.5121273941235152		-0.6088862173867386
 0.6484150302687719		-0.5255024588725531
 4.561668249671571		-0.39108978096909386
 5.567667636297536		-0.47736482122261414]';

Q_12 = 2*atan(x_12);

kin_double.joint_type = kin.joint_type;
kin_double.H = double(kin.H);
kin_double.P = double(kin.P);

Q = polynomial_IK.general_6R_given_q12(kin_double, R_06, T, Q_12)

%% Test solutions
q = Q(:,8)
[R_t, p_t] = fwdkin(kin_double, q)
R_t - R_06
p_t - T

%% Make gif showing all solutions

filename = "yumi_solns.gif";
clear im
for i = 1:8
h_fig = diagrams.setup;
hold on
options = {"show_arrows", true,...
           "unit_size", 40,...
           "show_arrow_labels", false,...
           "show_joint_labels", false,...
           "show_base_label", false,...
           "show_base_frame", false,...
           "show_task_frame", false,...
           "cyl_half_length", 20,...
           "cyl_radius", 10};
diagrams.robot_plot(kin_double, Q(:,i), options{:});
campos([-2e3 -2e3 2e3])
camva(8)
camtarget([-50 -250 200])

hold off
diagrams.redraw;
frame = getframe(h_fig);
im{i} = frame2im(frame);
end

for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.5);
    else
        imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.5);
    end
end
