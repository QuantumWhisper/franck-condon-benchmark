function result = mFCMatrix(varargin)

temp = memoize(@FCMatrix);
temp.CacheSize = 1e+6;
result = temp(varargin{:});

end