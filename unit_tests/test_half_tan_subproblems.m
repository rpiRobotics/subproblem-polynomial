%%
N_tests = 1
%%
N_tests = 1e4
%% SP1 direct to R
setup = subproblem_setups.sp_1;
e_vec = NaN(1, N_tests);
for i = 1:N_tests
    [P, S_given] = setup.setup();

    R = half_tan_sp1_R(P.p1, P.p2, P.k);
    e_vec(i) = norm(R*P.p1 - P.p2);
end


%% SP1 comparison
setup = subproblem_setups.sp_1;
e_vec = NaN(1, N_tests);
theta_vec = NaN(1, N_tests);
for i = 1:N_tests
    [P, S_given] = setup.setup();
    S.theta = subproblem.sp_1(P.p1, P.p2, P.k);
    e_vec(i) = setup.error(P,S);
    theta_vec(i) = S.theta;
end
%% SP1

setup = subproblem_setups.sp_1;
e_vec = NaN(1, N_tests);
theta_vec = NaN(1, N_tests);
x_vec = NaN(1, N_tests);
for i = 1:N_tests
    [P, S_given] = setup.setup();

    x = half_tan_sp1(P.p1, P.p2, P.k);
    S.theta = 2*atan(x);
    e_vec(i) = setup.error(P,S);
    theta_vec(i) = S.theta;
    x_vec(i) = x;
end

%% SP2
setup = subproblem_setups.sp_2;
e_vec = NaN(1, N_tests);
for i = 1:N_tests
    [P, S_given] = setup.setup();

    [x1, x2, a1, b1, c1, a2, b2, c2] = half_tan_sp2(P.p1, P.p2, P.k1, P.k2);
    S.theta1 = 2*atan(x1);
    S.theta2 = 2*atan(x2);
    e_vec(i) = setup.error(P,S);
end

%% SP3
setup = subproblem_setups.sp_3;
e_vec = NaN(1, N_tests);
for i = 1:N_tests
    [P, S_given] = setup.setup();

    [x, a, b, c] = half_tan_sp3(P.p1, P.p2, P.k, P.d);
    S.theta = 2*atan(x);
    e_vec(i) = setup.error(P,S);
end

%% SP4
setup = subproblem_setups.sp_4;
e_vec = NaN(1, N_tests);
for i = 1:N_tests
    [P, S_given] = setup.setup();

    [x, a, b, c] = half_tan_sp4(P.h, P.p, P.k, P.d);
    S.theta = 2*atan(x);
    e_vec(i) = setup.error(P,S);
end

%% Display results
semilogy(sort(e_vec))
yline(mean(e_vec));

%%
semilogy(theta_vec, e_vec, '.')
xlabel("\theta")
ylabel("Subproblem 1 Error")

%%
loglog(abs(x_vec), e_vec, '.')
xlabel("|x|")
ylabel("Subproblem 1 Error")