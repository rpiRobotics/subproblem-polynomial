function Q = three_pairs_intersecting_given_q4(kin, R_06, p_0T, Q_4)
Q = [];
is_LS_vec = [];
for i_q4 = 1:length(Q_4)
    [e, Q_partial] = alignment_err_given_q4(Q_4(i_q4), p_0T, R_06, kin);
    soln_num = abs(e) < 1e-6;
    q_partial = Q_partial(:,soln_num);

    R_04 = rot(kin.H(:,1),q_partial(1)) * rot(kin.H(:,2),q_partial(2)) ...
         * rot(kin.H(:,3),q_partial(3)) * rot(kin.H(:,4),q_partial(4));
    [q5, q5_is_LS] = subproblem.sp_1(kin.H(:,6), R_04'*R_06*kin.H(:,6),  kin.H(:,5));
    [q6, q6_is_LS] = subproblem.sp_1(kin.H(:,5), R_06'*R_04*kin.H(:,5), -kin.H(:,6));
    q_i = [q_partial; q5; q6];
    Q = [Q q_i];
    is_LS_vec = [is_LS_vec [q5_is_LS; q6_is_LS; e(soln_num)] ];
end
end

function [e_vec, Q_partial] = alignment_err_given_q4(q4, p_16, R_06, kin)
    e_vec = NaN([1 4]);
    Q_partial = NaN([4 4]);
    i_soln = 1;
    % Find 2 solutions of q3 with Subproblem 3
    % Find 2 solutions of (q1, q2) with Subproblem 2
    R_34 = rot(kin.H(:,4), q4);
    
    [t3, t3_is_LS] = subproblem.sp_3(R_34*kin.P(:,5), -kin.P(:,3), kin.H(:,3), norm(p_16));
    if t3_is_LS
        return
    end

    for i_t3 = 1:length(t3)
        q3 = t3(i_t3);
        R_23 = rot(kin.H(:,3), q3);

        [t1, t2, t12_is_LS] = subproblem.sp_2(p_16, kin.P(:,3)+R_23*R_34*kin.P(:,5), -kin.H(:,1), kin.H(:,2));
        if t12_is_LS
            i_soln = i_soln + 2;
            continue
        end
        for i_12 = 1:length(t1)
            q1 = t1(i_12);
            q2 = t2(i_12);
            R_01 = rot(kin.H(:,1), q1);
            R_12 = rot(kin.H(:,2), q2);
            R_04 = R_01 * R_12 * R_23 * R_34;
            e_i = kin.H(:,5)' * R_04' * R_06* kin.H(:,6) - kin.H(:,5)'*kin.H(:,6);
            e_vec(i_soln) = e_i;
            Q_partial(:,i_soln) = [q1; q2; q3; q4];
            i_soln = i_soln + 1;
        end
    end
end