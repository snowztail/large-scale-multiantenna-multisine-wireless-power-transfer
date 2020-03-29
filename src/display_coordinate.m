function [coordinateString] = display_coordinate(coordinateMatrix)
    % Function:
    %   - convert the coordinate from matrix form to conventional notation
    %
    % InputArg(s):
    %   - coordinateMatrix: rows denote axes, columns denote entries
    %
    % OutputArg(s):
    %   - coordinateString: of human-readable format (x,y) or (x,y,z)
    %
    % Comment(s):
    %   - for 2-D and 3-D coordinates only
    %
    % Author & Date: Yang (i@snowztail.com) - 29 Mar 20


    nCoordinates = size(coordinateMatrix, 2);
    coordinateString = cell(nCoordinates, 1);
    for iCoordinate = 1 : nCoordinates
        if size(coordinateMatrix, 1) == 2
            coordinateString{iCoordinate} = sprintf('(%d,%d)', coordinateMatrix(1, iCoordinate), coordinateMatrix(2, iCoordinate));
        elseif size(coordinateMatrix, 1) == 3
            coordinateString{iCoordinate} = sprintf('(%d,%d,%d)', coordinateMatrix(1, iCoordinate), coordinateMatrix(2, iCoordinate), coordinateMatrix(3, iCoordinate));
        else
            error('Sorry, this function currently supports 2-D and 3-D inputs only...');
        end
    end
end
