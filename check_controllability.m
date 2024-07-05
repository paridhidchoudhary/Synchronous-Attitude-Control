function is_controllable = check_controllability(A, B)
    % Check the controllability of the system
    Co = ctrb(A, B);
    rank_Co = rank(Co);
    is_controllable = (rank_Co == size(A, 1));
end