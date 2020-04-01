function [coefficient] = rr_equations(delta, component, constraint, nTxs, nSubbands, nUsers, userIndex)
    % Function:
    %   - form the equations that the solution delta should satisfy
    %
    % InputArg(s):
    %   - delta [\Delta] (waveformRank * waveformRank): the Hermitian matrix used in rank reduction
    %   - component [V] ((nTxs * nSubbands) * waveformRank): component of the waveform matrix such that V times its Hermitian equals X
    %   - constraint [A] (1, nUsers) {(nTxs * nSubbands) * (nTxs * nSubbands)}: QCQP constraints
    %   - nTxs [M]: number of transmit antennas
    %   - nSubbands [N]: number of subbands/subcarriers
    %   - nUsers [K]: number of users
    %   - userIndex [q_0]: index of user with the largest output voltage
    %
    % OutputArg(s):
    %   - coefficient: the coefficients corresponding to equations as a function of delta
    %
    % Comment(s):
    %   - all the equations are functions of delta and should equal 0
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %   - Y. Huang and D. Palomar, "Rank-Constrained Separable Semidefinite Programming With Applications to Optimal Beamforming," IEEE Transactions on Signal Processing, vol. 58, no. 2, pp. 664–678, 2010.
    %
    % Author & Date: Yang (i@snowztail.com) - 1 Apr 20



    % * the first set is target constraint: the target functions Tr{\boldsymbol{AX}} + \bar{C} of all users are no greater than that of the indexed user
    targetCoefficient = zeros(nUsers, 1);
    for iUser = 1 : nUsers
        if iUser ~= userIndex
            targetCoefficient(iUser) = real(trace(component' * (constraint{iUser} - constraint{userIndex}) * component * delta));
        end
    end
    targetCoefficient(userIndex) = [];

    % * the second set is transmit power constraint: the power of waveform Tr{X} is no greater than budget
    powerCoefficient = real(trace(component' * eye(nTxs * nSubbands) * component * delta));

    % * the third set is Hermitian constraint: delta is Hermitian matrix
    hermitianCoefficient = real(vec(delta') - vec(delta));

    % * combine coefficients
    coefficient = [targetCoefficient; powerCoefficient; hermitianCoefficient];

end
