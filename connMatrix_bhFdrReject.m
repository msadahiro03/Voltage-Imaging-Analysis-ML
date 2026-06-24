function rejectMask = connMatrix_bhFdrReject(pValues, q)
%CONNMATRIX_BHFDRREJECT Benjamini-Hochberg FDR rejection mask.
%   rejectMask = connMatrix_bhFdrReject(pValues, q)
%   pValues: any shape; NaN/non-finite entries are excluded from the FDR family
%   and are never rejected. q: FDR level in (0, 1].

rejectMask = false(size(pValues));
fin = isfinite(pValues);
if ~any(fin(:))
    return
end

if ~(isscalar(q) && isnumeric(q) && isfinite(q) && q > 0 && q <= 1)
    error('connMatrix_bhFdrReject: q must be a scalar in (0, 1].');
end

p = pValues(fin);
m = numel(p);
if m == 0
    return
end

[sortedP, sortIdx] = sort(p);
k = 0;
for ii = m:-1:1
    if sortedP(ii) <= (ii / m) * q
        k = ii;
        break
    end
end

if k == 0
    return
end

r = false(m, 1);
r(sortIdx(1:k)) = true;
rejectMask(fin) = r;
end
