function M = FCMatrix(q1,q2,lambda)

nq2 = numel(q2);
nq1 = numel(q1);
nl = numel(lambda);
M = zeros([nq1,nq2,nl]);
for ll = 1:nl
    for ii = 1:nq1
        M(ii,:,ll) = FCMatrixMulti(q1(ii),q2,lambda(ll));
    end
end
if nq1 == 1 && nq2 == 1
    M = M(:);
end
if iscolumn(M)
    M = M';
end


end

function M = FCMatrixMulti(q1,q2,lambda)

nq2 = numel(q2);
nq1 = numel(q1);
nl = numel(lambda);

if nq1 > 1
    M = zeros(size(q1));
elseif nq2 > 1
    M = zeros(size(q2));
elseif nl > 1
    M = zeros(size(lambda));
end

M(:) = mFCMatrixSingle(q1,q2,lambda);

end

function result = mFCMatrixSingle(varargin)

temp = memoize(@FCMatrixSingle);
temp.CacheSize = 1e+8;
result = temp(varargin{:});

end
