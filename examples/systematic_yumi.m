% kin = hardcoded_IK_setups.yumi_fixed_q3.get_kin
%%
% pi_sym = sym(pi);
% 
% alpha_vec = 2*pi_sym/360*[-90 90 -90 -90 -90 90 0];
% a_vec     = [-30 30 40.5 40.5 27 -27 0];
% d_vec     = [166 0 251.5 0 265 0 36];
% 
% kin_full = dh_to_kin(alpha_vec, a_vec, d_vec);
% 
% [kin, R_6T] = fwdkin_partial(kin_full, pi_sym/2, 3);
% 
% kin_full.P(:,end) = 0;
% 
% [R_t, p_t] = fwdkin(kin_full, 2*pi_sym/360*([0 0 0 -90 180 0 0]))
%%
pi_sym = sym(pi);
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
zv = [0;0;0];

kin_full.joint_type = zeros([1 7]);
kin_full.H = [ez ey ez ey ex ey ex];
kin_full.P = [zv 166*ez-30*ex 30*ex 251.5*ez+40.5*ex 40.5*ez 265*ex-27*ez 27*ez zv];

% [R_t, p_t] = fwdkin(kin_full, zeros([7 1]))

[kin, R_6T] = fwdkin_partial(kin_full, pi_sym/2, 3);

%% Instead, fix q1

[kin, R_6T] = fwdkin_partial(kin_full, 0, 1);
kin.P(:,1) = 0;

%% Different q_3

[kin, R_6T] = fwdkin_partial(kin_full, pi_sym/6, 3);
kin.P = rot(ez,-pi_sym/6)*kin.P;
kin.H = rot(ez,-pi_sym/6)*kin.H;
%%
% R_06 = rot([0 1 0]', pi_sym/2)
% T = [200 200 200]'
%% Plot robot pose
diagrams.setup;
hold on
options = {"show_arrows", true,...
           "unit_size", 20,...
           "show_arrow_labels", false,...
           "show_joint_labels", false,...
           "show_base_label", false,...
           "show_base_frame", false,...
           "show_task_frame", false,...
           "cyl_half_length", 10,...
           "cyl_radius", 5};

% diagrams.robot_plot(kin_full, deg2rad([0 0 0 -90 180 0 0]), options{:});
% diagrams.robot_plot(kin_full, zeros([7 1]), options{:});
% diagrams.robot_plot(kin, deg2rad([0 0 -90 180 0 0]), options{:});
diagrams.robot_plot(kin_double, zeros([6 1]), options{:});

hold off
diagrams.redraw;
%% Make sure search-based IK works
close
R_06 = eye(3);
% T = [300; 70; 160];
T = [200 200 200]';
% T = rand_vec*500;
% T =[  254
%  -130
%    68]

% T =[50
%  -355
%   354]

% [R_06, T] = fwdkin(kin_double, rand_angle([6 1]))

kin_double.joint_type = kin.joint_type;
kin_double.H = double(kin.H);
kin_double.P = double(kin.P);

[Q, is_LS_vec] = IK.IK_gen_6_dof_mex(double(R_06), double(T), kin_double)


[R_06_t, T_t] = fwdkin(kin_double, Q(:,1))
view(0,-90)
%%
p_0T = T;
p_06 = p_0T - kin.P(:,1) - R_06*kin.P(:,7);

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
% R_05 = R_01 * R_12 * R_23 * R_34 * R_45;
%lhs_err = R_05 * kin.H(:,6) - R_06 * kin.H(:,6); % 3 equations
%lhs_err = (R_06 * kin.H(:,6))'*  R_05 * kin.H(:,6) -1
lhs_err_1 = kin.H(:,3)'*R_34*R_45*kin.H(:,6) - kin.H(:,3)'*(R_01*R_12)'*R_06*kin.H(:,6); % No x3
lhs_err_2 = kin.H(:,4)'*R_45*kin.H(:,6) - kin.H(:,4)'*(R_01*R_12*R_23)'*R_06*kin.H(:,6); % No x4

% lhs_err_3 = kin.H(:,1)'*R_12 * R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,1)'*R_06*kin.H(:,6); % no x1
% lhs_err_4 = kin.H(:,2)'*R_23 * R_34 * R_45 * kin.H(:,6) - kin.H(:,2)'* R_01' * R_06*kin.H(:,6); % no x2
% lhs_err_5 = kin.H(:,5)'*kin.H(:,6) - kin.H(:,5)'*(R_01*R_12*R_23*R_34)'*R_06*kin.H(:,6); % No x5


eqns_lhs = [lhs_err_1; lhs_err_2; lhs_1; lhs_2];

eqns_lhs = simplify(eqns_lhs)


eqns_frac = simplifyFraction(eqns_lhs)
eqns_num = numden(eqns_frac)
%%
fileID = fopen('yumi_eqns.txt','w');
for i = 1:length(eqns_num)
    fprintf(fileID, "p%i = ", i);
    % fprintf(fileID, string(vpa(eqns_num(i), 6)) + '\n');
    fprintf(fileID, string(eqns_num(i)) + '\n');
end

%% Test eqns

q = Q(:,1);
x = tan(q/2)

double(subs(eqns_lhs(4), [x1 x2 x3 x5], x([1 2 3 5])'))
double(subs(eqns_num(4), [x1 x2 x3 x5], x([1 2 3 5])'))

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

Q = [];
for i = 1:8
q_12 = Q_12(:,i)
[e_vec, Q_i, is_LS_vec] = alignment_err_given_q12(q_12(1), q_12(2), T, R_06, kin_double)
i_soln = e_vec<1e-6;
Q = [Q Q_i(:,i_soln)];
end
disp(Q)

%%
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
%%
q = Q(:,8)
[R_t, p_t] = fwdkin(kin_double, q)
R_t - R_06
p_t - T

function [e_vec, Q, is_LS_vec] = alignment_err_given_q12(q1, q2, p_16, R_06, kin)
    e_vec = NaN([1 4]);
    Q = NaN([6 4]);
    is_LS_vec = NaN([2 4]);
    % find up to 4 solutions of (q3, q4, q5) using Subproblem 5
    p_63 = rot(-kin.H(:,2),q2) * (rot(-kin.H(:,1),q1)*p_16 - kin.P(:,2)) - kin.P(:,3);
    [t3, t4, t5] = subproblem.sp_5( ...
                   -kin.P(:,4), p_63, kin.P(:,5), kin.P(:,6), ...
                   -kin.H(:,3), kin.H(:,4), kin.H(:,5));
    for i_q345 = 1:length(t3)
        q3 = t3(i_q345);
        q4 = t4(i_q345);
        q5 = t5(i_q345);

        % Calculate alignment error for h_6
        R_05 = rot(kin.H(:,1),q1)*rot(kin.H(:,2),q2) ...
            * rot(kin.H(:,3),q3)*rot(kin.H(:,4),q4) ...
            * rot(kin.H(:,5),q5);

        e_vec(i_q345) = norm(R_05*kin.H(:,6)-R_06*kin.H(:,6));

        if nargout > 1
            % Find q6 with subproblem 1
            p = kin.H(:,5); % Can't be collinear with h_6
            [q6, q6_is_ls] = subproblem.sp_1(p, R_05'*R_06*p, kin.H(:,6));
            % TODO: q6_is_ls will tend to be true since e_i isn't exactly 0

            Q(:,i_q345) = [q1; q2; q3; q4; q5; q6];
            is_LS_vec(:,i_q345) = [q6_is_ls; e_vec(i_q345)];
        end
    end
end