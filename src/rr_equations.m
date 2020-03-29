function [value] = rr_equations(delta, component, term, nTxs, nSubbands, nUsers, userIndex)



    % delta = cell(nUsers, 1);

    % for iUser = 1 : nUsers
    %     delta{iUser} = reshape(deltaVector(1 + sum(precoderRank(1 : iUser - 1) .^ 2) : sum(precoderRank(1 : iUser) .^ 2)), precoderRank(iUser), precoderRank(iUser));
    % end

    % trace constraints on upperbound
    value1 = zeros(nUsers, 1);
    for iUser = 1 : nUsers
        if iUser ~= userIndex
            value1(iUser) = real(trace(component' * (term{iUser} - term{userIndex}) * component* delta));
        end
    end
    value1(userIndex) = [];

    % trace constraints on transmit power
    value2 = real(trace(component' * eye(nTxs * nSubbands) * component * delta));

    % Hermitian constraints
    value3 = real(vec(delta') - vec(delta));

    value = [value1; value2; value3];

end
