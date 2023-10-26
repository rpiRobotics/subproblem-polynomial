function Q = general_6R_given_q12(kin, R_06, p_06, Q_12)
Q = [];
for i = 1:width(Q_12)
    q_12 = Q_12(:,i);
    [e_vec, Q_i, is_LS_vec] = alignment_err_given_q12(q_12(1), q_12(2), p_06, R_06, kin);
    i_soln = e_vec<1e-6;
    Q = [Q Q_i(:,i_soln)];
end
end

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