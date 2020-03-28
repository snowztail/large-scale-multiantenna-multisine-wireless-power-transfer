function [value] = equations(deltaVector, factor, term, pathloss, nSubbands, nUsers, userIndex, precoderRank)



    delta = cell(nUsers, 1);

    for iUser = 1 : nUsers
        delta{iUser} = reshape(deltaVector(1 + sum(precoderRank(1 : iUser - 1) .^ 2) : sum(precoderRank(1 : iUser) .^ 2)), precoderRank(iUser), precoderRank(iUser));
    end

    % trace constraints on upperbound
    value1 = zeros(nUsers, 1);
    for iUser = 1 : nUsers
        if iUser ~= userIndex
            value1(iUser) = real(trace(factor{iUser}' * term{iUser} * factor{iUser} * delta{iUser})) - real(trace(factor{userIndex}' * term{userIndex} * factor{userIndex} * delta{userIndex}));
        end
    end
    value1(userIndex) = [];

    % trace constraints on transmit power
    value2 = 0;
    for iUser = 1 : nUsers
        value2 = value2 + real(trace(factor{iUser}' * (pathloss(iUser) * eye(nSubbands)) * factor{iUser} * delta{iUser}));
    end

    value3 = zeros(length(deltaVector), 1);
    % Hermitian constraints
    for iUser = 1 : nUsers
        value3(1 + sum(precoderRank(1 : iUser - 1) .^ 2) : sum(precoderRank(1 : iUser) .^ 2)) = real(vec(delta{iUser}') - vec(delta{iUser}));
    end

    value = [value1; value2; value3];

end
