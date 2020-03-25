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
    %   - assume the input matrix is Hermitian
    %   - the output component is not unique
    %
    % Reference(s):
    %   - G. H. Golub and C. F. van. Loan, Matrix computations. Baltimore: Johns Hopkins University Press, 2013.
    %
    % Author & Date: Yang (i@snowztail.com) - 25 Mar 20


    l = length(matrix);
    r = rank(matrix);
    [v, d] = eig(matrix);
    [d, index] = sort(diag(d), 'ascend', 'ComparisonMethod', 'abs');
    v = v(:, index);
    component = v * [zeros(l - r, r); sqrt(diag(d(l - r + 1 : end)))];

end
