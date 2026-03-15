function w = rateW_lead(rateW4DCache,N,n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,Vg)
% Computes total rate summed over both leads (L and R).
% Used by generateMatrixW for the rate equation matrix.

lead = 1;
temp_wl = rateW(rateW4DCache,n1,n2,q1,0:(N-1),vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg);
wl = temp_wl(q2+1);
lead = -1;
temp_wr = rateW(rateW4DCache,n1,n2,q1,0:(N-1),vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg);
wr = temp_wr(q2+1);
w = wl + wr;

end
