syms x real
syms D real
R_67 = half_tan_rot(kin.H(:,7), x);

p_07 = p_0T - R_07 * kin.P(:,8) + [0;0;D];
[~, n_SEW] = SEW.inv_kin([0;0;0], p_07, psi);

% q6 (SP4 based on q7 search)
[x6, a, b, c] = half_tan_sp4(R_67*R_07'*n_SEW, kin.P(:,6), -kin.H(:,6), n_SEW'*p_07);
x_6 = simplify(x6(1))

%%
syms x_6 real
R_56 = half_tan_rot(kin.H(:,6), x_6);
p_WE = R_07 * R_67' * R_56' * kin.P(:,6);

% q1
[x1, a, b, c] = half_tan_sp3(kin.P(:,2), p_07-p_WE, kin.H(:,1), norm(kin.P(:,4)));
x_1 = simplify(x1(1))
%%
syms x_1 real
R_01 = half_tan_rot(kin.H(:,1), x_1);

[x3, x2, a3, b3, c3, a2, b2, c2] = half_tan_sp2(kin.P(:,4), R_01'*p_07-R_01'*p_WE-kin.P(:,2), kin.H(:,3), -kin.H(:,2));
x2 = simplify(x2(1))
x3 = simplify(x3(1))
%%
syms x_2 x_3 real
R_12 = half_tan_rot(kin.H(:,2), x_2);
R_23 = half_tan_rot(kin.H(:,3), x_3);

% Error eqn
lhs_err = kin.H(:,4)' * (R_01 * R_12 * R_23)' * R_07 * (R_56 * R_67)'* kin.H(:,5) - kin.H(:,4)' * kin.H(:,5);
lhs_err = simplify(lhs_err)

%%
%% Plug into DESMOS
x_6 = -(8000*x + (14142231*x^4 + 99715538*x^2 + 14142231)^{1/2}(S_1))/(2637*x^2 + 5363)

x_1 = -((77*(x_6^2 - 1))/(500*(x_6^2 + 1)) - (110403*(x^2 - 1))/(5000000*(x^2 + 1)) + (((77*(x_6^2 - 1))/(500*(x_6^2 + 1)) - (110403*(x^2 - 1))/(5000000*(x^2 + 1)) + (162*x*x_6)/(625*(x^2 + 1)*(x_6^2 + 1)) + 223/2000)^2 - ((104951*(x^2 - 1))/(1000000*(x^2 + 1)) + (81*(x_6^2 - 1))/(625*(x_6^2 + 1)) + 2*((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/2)^2 + 2*(D + (1363*x)/(5000*(x^2 + 1)) + (4*x_6*(x^2 - 1))/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/4)^2 + 2*((2*(x_6^2 - 1))/(5*(x_6^2 + 1)) + 1/2)^2 - (154*x*x_6)/(125*(x^2 + 1)*(x_6^2 + 1)) + 289862931016633453/1125899906842624000)*(((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/2)^2/2 - (81*(x_6^2 - 1))/(2500*(x_6^2 + 1)) - (104951*(x^2 - 1))/(4000000*(x^2 + 1)) + (D + (1363*x)/(5000*(x^2 + 1)) + (4*x_6*(x^2 - 1))/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/4)^2/2 + ((2*(x_6^2 - 1))/(5*(x_6^2 + 1)) + 1/2)^2/2 + (77*x*x_6)/(250*(x^2 + 1)*(x_6^2 + 1)) - 941871567069197203/4503599627370496000))^{1/2}(S_2) + (162*x*x_6)/(625*(x^2 + 1)*(x_6^2 + 1)) + 223/2000)/((104951*(x^2 - 1))/(2000000*(x^2 + 1)) + (81*(x_6^2 - 1))/(1250*(x_6^2 + 1)) + ((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/2)^2 + (D + (1363*x)/(5000*(x^2 + 1)) + (4*x_6*(x^2 - 1))/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/4)^2 + ((2*(x_6^2 - 1))/(5*(x_6^2 + 1)) + 1/2)^2 - (77*x*x_6)/(125*(x^2 + 1)*(x_6^2 + 1)) + 289862931016633453/2251799813685248000)
 
x_2 = -(2*D - ((2*D + (1363*x)/(2500*(x^2 + 1)) + (8*x_6*(x^2 - 1))/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/2)^2 - ((4*x_1)/(x_1^2 + 1) - (2*(x_1^2 - 1))/(x_1^2 + 1) + (8*x_1*((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1))))/(x_1^2 + 1) - (8*(x_1^2 - 1)*(x_6^2 - 1))/(5*(x_1^2 + 1)*(x_6^2 + 1)) + 319/250)*((x_1^2 - 1)/(2*(x_1^2 + 1)) - x_1/(x_1^2 + 1) - (2*x_1*((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1))))/(x_1^2 + 1) + (2*(x_1^2 - 1)*(x_6^2 - 1))/(5*(x_1^2 + 1)*(x_6^2 + 1)) + 481/1000))^{1/2}(S_3) + (1363*x)/(2500*(x^2 + 1)) + (8*x_6*(x^2 - 1))/(5*(x^2 + 1)*(x_6^2 + 1)) + 1/2)/((2*x_1)/(x_1^2 + 1) - (x_1^2 - 1)/(x_1^2 + 1) + (4*x_1*((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1))))/(x_1^2 + 1) - (4*(x_1^2 - 1)*(x_6^2 - 1))/(5*(x_1^2 + 1)*(x_6^2 + 1)) + 319/500)
  
x_3 = -(-((x_1^2 - 1)/(2*(x_1^2 + 1)) + x_1/(x_1^2 + 1) + (((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1)))*(x_1^2 - 1))/(x_1^2 + 1) + (4*x_1*(x_6^2 - 1))/(5*(x_1^2 + 1)*(x_6^2 + 1)) + 3/125)*((2*(x_1^2 - 1))/(x_1^2 + 1) + (4*x_1)/(x_1^2 + 1) + (4*((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1)))*(x_1^2 - 1))/(x_1^2 + 1) + (16*x_1*(x_6^2 - 1))/(5*(x_1^2 + 1)*(x_6^2 + 1)) + 361/250))^{1/2}(S_4)/((x_1^2 - 1)/(x_1^2 + 1) + (2*x_1)/(x_1^2 + 1) + (2*((1363*(x^2 - 1))/(10000*(x^2 + 1)) - (8*x*x_6)/(5*(x^2 + 1)*(x_6^2 + 1)))*(x_1^2 - 1))/(x_1^2 + 1) + (8*x_1*(x_6^2 - 1))/(5*(x_1^2 + 1)*(x_6^2 + 1)) + 361/500)
 
e =  (4*x*x_6*(((x_1^2 - 1)*(x_3^2 - 1))/((x_1^2 + 1)*(x_3^2 + 1)) + (8*x_1*x_2*x_3)/((x_1^2 + 1)*(x_2^2 + 1)*(x_3^2 + 1))))/((x^2 + 1)*(x_6^2 + 1)) - ((x_6^2 - 1)*((2*x_1*(x_3^2 - 1))/((x_1^2 + 1)*(x_3^2 + 1)) - (4*x_2*x_3*(x_1^2 - 1))/((x_1^2 + 1)*(x_2^2 + 1)*(x_3^2 + 1))))/(x_6^2 + 1) + (4*x_3*x_6*(x^2 - 1)*(x_2^2 - 1))/((x^2 + 1)*(x_2^2 + 1)*(x_3^2 + 1)*(x_6^2 + 1))
 