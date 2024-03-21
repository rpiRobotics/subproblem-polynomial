function [Q, is_LS_vec] = sawyer_find_remaining_angles(R_07, p_0T, SEW_class, psi, kin, q7_vec)
% Almost same as IK_R_2R_2R_2R_q7, but 1D search is already solved
% Just need to check which branch is correct
THRESH = 1e-3;

Q = [];
is_LS_vec = [];

for i = 1:length(q7_vec)
    [alignment, q_solns_partial_i] = alignment_given_q7(q7_vec(i), kin, R_07, p_0T, psi, SEW_class);
    % q_partial_i = q_solns_partial_i(:,soln_num_vec(i));
    for j = 1:width(q_solns_partial_i)
        if isnan(alignment(j)) || abs(alignment(j)) >= THRESH
            continue
        end
        [q, is_LS] = q_given_q123__67(q_solns_partial_i(:,j), kin, R_07);
        Q = [Q q];
        is_LS_vec = [is_LS_vec is_LS];
    end
end
end


function [alignment, q_solns_partial] = alignment_given_q7(q7, kin, R_07, p_0T, psi, SEW)
    q_solns_partial = NaN(7, 8);
    alignment = NaN(1, 8);
    i_soln = 1;

    p_07 = p_0T - R_07*kin.P(:,8);
    R_67 = rot(kin.H(:,7), q7);
    [~, n_SEW] = SEW.inv_kin([0;0;0], p_07, psi);

    % Find q_6 with Subproblem 4
    [t6, t6_is_LS] = subproblem.sp_4(R_67*R_07'*n_SEW, kin.P(:,6), -kin.H(:,6), n_SEW'*p_07);
    if t6_is_LS
        return
    end

    for i_t6 = 1:length(t6)
        q6 = t6(i_t6);
        R_56 = rot(kin.H(:,6), q6);
        p_WE = R_07 * R_67' * R_56' * kin.P(:,6);

        % Find q_1 with Subproblem 3
        [t1, t1_is_LS] = subproblem.sp_3(kin.P(:,2), p_07-p_WE, kin.H(:,1), norm(kin.P(:,4)));
        if t1_is_LS
            i_soln = i_soln + 4; % TODO check 
            continue
        end

        for i_t1 = 1:length(t1)
            q1 = t1(i_t1);
            R_01 = rot(kin.H(:,1), q1);

            % Find (q_2, q_3) with Subproblem 2
            [t3, t2, t23_is_LS] = subproblem.sp_2(kin.P(:,4), R_01'*p_07-R_01'*p_WE-kin.P(:,2), kin.H(:,3), -kin.H(:,2));
            if t23_is_LS
                i_soln = i_soln + 2; % TODO check 
                continue
            end

            for i_23 = 1:length(t2)
                q2 = t2(i_23);
                q3 = t3(i_23);
                R_12 = rot(kin.H(:,2), q2);
                R_23 = rot(kin.H(:,3), q3);
                R_03 = R_01*R_12*R_23;
                R_57 = R_56 * R_67;

                e_i = kin.H(:,4)' * R_03' * R_07 * R_57' * kin.H(:,5) - kin.H(:,4)'* kin.H(:,5);
                alignment(i_soln) = e_i;
                q_solns_partial(:,i_soln) = [q1; q2; q3; NaN; NaN; q6; q7];
                i_soln = i_soln + 1;
            end
        end
    end
end


function [q, is_LS] = q_given_q123__67(q123__67, kin, R_07)    
    R_01 = rot(kin.H(:,1), q123__67(1));
    R_12 = rot(kin.H(:,2), q123__67(2));
    R_23 = rot(kin.H(:,3), q123__67(3));
    R_56 = rot(kin.H(:,6), q123__67(6));
    R_67 = rot(kin.H(:,7), q123__67(7));

    R_03 = R_01*R_12*R_23;
    R_57 = R_56*R_67;

    [q4, q4_is_LS] = subproblem.sp_1(kin.H(:,5), R_03'*R_07 *R_57'*kin.H(:,5), kin.H(:,4));
    [q5, q5_is_LS] = subproblem.sp_1(kin.H(:,4), R_57 *R_07'*R_03 *kin.H(:,4), -kin.H(:,5));

    q = [q123__67(1:3); q4; q5; q123__67(6:7)];
    is_LS = q4_is_LS || q5_is_LS;
end
