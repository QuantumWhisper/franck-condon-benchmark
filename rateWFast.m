function w = rateWFast(rateW4D,n1,n2,q1,q2,varargin)

w = rateW4D(q1+1,q2+1,n1+1,n2+1);
w(isnan(w)) = 0;

end