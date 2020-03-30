function [component] = decompose(matrix)
    % Function:
    %   - decompose a Hermitian matrix into the product of a component and its Hermitian
    %
    % InputArg(s):
    %   - matrix (M * M): a Hermitian matrix with rank R
    %
    % OutputArg(s):
    %   - component (M * R): a candidate component
    %
    % Comment(s):
    %   - for Hermitian inputs only
    %   - the output component is not unique
    %
    % Reference(s):
    %   - G. H. Golub and C. F. van. Loan, Matrix computations. Baltimore: Johns Hopkins University Press, 2013.
    %
    % Author & Date: Yang (i@snowztail.com) - 25 Mar 20


    r = rank(matrix);
    [u, s] = svd(matrix);
    [s, index] = sort(diag(s), 'descend', 'ComparisonMethod', 'abs');
    u = u(:, index);
    s = diag(s);
    component = u(:, 1 : r) * s(1 : r, 1 : r) ^ (1 / 2);

end
