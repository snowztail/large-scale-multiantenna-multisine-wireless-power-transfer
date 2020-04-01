function [coefficient] = rr_equations_che(deltaVector, component, constraint, pathloss, nSubbands, nUsers, userIndex, carrierWeightRank)
    % Function:
    %   - form the equations that the solution delta should satisfy
    %
    % InputArg(s):
    %   - deltaVector [\Delta] (sum(carrierWeightRank .^ 2) * 1): vectorized form of the Hermitian matrix used in rank reduction
    %   - component [V] (nUsers, 1) {(nTxs * nSubbands) * carrierWeightRank}: component of the carrier weight matrices such that V times its Hermitian equals X
    %   - constraint [A] (1, nUsers) {(nTxs * nSubbands) * (nTxs * nSubbands)}: QCQP constraints
    %   - pathloss [\boldsymbol{\Lambda}] (1 * nUsers) : large-scale channel strength reduction
    %   - nSubbands [N]: number of subbands/subcarriers
    %   - nUsers [K]: number of users
    %   - userIndex [q_0]: index of user with the largest output voltage
    %   - carrierWeightRank [R] (nUsers, 1): rank of carrier weight matrices
    %
    % OutputArg(s):
    %   - coefficient: the coefficients corresponding to equations as a function of delta
    %
    % Comment(s):
    %   - all the equations are functions of delta and should equal 0
    %   - fsolve cannot handle cells and we reform the deltas as a column vector
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %   - Y. Huang and D. Palomar, "Rank-Constrained Separable Semidefinite Programming With Applications to Optimal Beamforming," IEEE Transactions on Signal Processing, vol. 58, no. 2, pp. 664–678, 2010.
    %
    % Author & Date: Yang (i@snowztail.com) - 1 Apr 20



    % * reform delta vector back to Hermitian matrices
    delta = cell(nUsers, 1);
    for iUser = 1 : nUsers
        delta{iUser} = reshape(deltaVector(1 + sum(carrierWeightRank(1 : iUser - 1) .^ 2) : sum(carrierWeightRank(1 : iUser) .^ 2)), carrierWeightRank(iUser), carrierWeightRank(iUser));
    end

    % * the first set is target constraint: the target functions Tr{\boldsymbol{AX}} + \bar{C} of all users are no greater than that of the indexed user
    targetCoefficient = zeros(nUsers, 1);
    for iUser = 1 : nUsers
        if iUser ~= userIndex
            targetCoefficient(iUser) = real(trace(component{iUser}' * constraint{iUser} * component{iUser} * delta{iUser})) - real(trace(component{userIndex}' * constraint{userIndex} * component{userIndex} * delta{userIndex}));
        end
    end
    targetCoefficient(userIndex) = [];

    % * the second set is transmit power constraint: the power of waveform Tr{X} is no greater than budget
    powerCoefficient = 0;
    for iUser = 1 : nUsers
        powerCoefficient = powerCoefficient + real(trace(component{iUser}' * (pathloss(iUser) * eye(nSubbands)) * component{iUser} * delta{iUser}));
    end

    % * the third set is Hermitian constraint: delta is Hermitian matrix
    hermitianCoefficient = zeros(length(deltaVector), 1);
    for iUser = 1 : nUsers
        hermitianCoefficient(1 + sum(carrierWeightRank(1 : iUser - 1) .^ 2) : sum(carrierWeightRank(1 : iUser) .^ 2)) = real(vec(delta{iUser}') - vec(delta{iUser}));
    end

    % * combine coefficients
    coefficient = [targetCoefficient; powerCoefficient; hermitianCoefficient];

end
