function w = m_rateW(n1,n2,q1,q2,vmode,alphaL,alphaR,lambda,Vsd,T,eta,lead,Vg)

gammaL = alphaL*vmode;
gammaR = alphaR*vmode;
epsilond = 0 + Vg;
muL = eta*Vsd;
muR = -(1-eta)*Vsd;
regularizedON = true;

if n1 == 1 && n2 == 0
    gamma = whichGamma(gammaL,gammaR,lead);
    mu = whichMu(muL,muR,lead);
    w = s_(n1,n2)*gamma/hbar_eV*(FCMatrix(q1,q2,lambda)).^2.*(1-fermi(epsilond-(q2-q1)*vmode,mu,T));
elseif n1 == 0 && n2 == 1
    gamma = whichGamma(gammaL,gammaR,lead);
    mu = whichMu(muL,muR,lead);
    w = s_(n1,n2)*gamma/hbar_eV*(FCMatrix(q1,q2,lambda)).^2.*(fermi(epsilond+(q2-q1)*vmode,mu,T));
elseif n1 == 0 && n2 == 0 
    if lead == 1% These transfer one electron from lead L to lead R,
        if regularizedON
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*(sumMMr(q1,q2,lambda,muL,muR,vmode,epsilond,T)+...
                sumMMMMrs(q1,q2,lambda,muL,muR,vmode,epsilond,T));
        else
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*integral(@(x) mmInt00(q1,q2,lambda,x,epsilond,vmode,T,muL,muR),-Inf,Inf);
        end
    else % These transfer one electron from lead R to lead L,
        if regularizedON
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*(sumMMr(q1,q2,lambda,muR,muL,vmode,epsilond,T)+...
                sumMMMMrs(q1,q2,lambda,muR,muL,vmode,epsilond,T));
        else
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*integral(@(x) mmInt00(q1,q2,lambda,x,epsilond,vmode,T,muR,muL),-Inf,Inf);
        end        
    end
elseif n1 == 1 && n2 == 1
    if lead == 1% These transfer one electron from lead L to lead R,
        if regularizedON
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*(sumMMr11(q1,q2,lambda,muL,muR,vmode,epsilond,T)+...
                sumMMMMrs11(q1,q2,lambda,muL,muR,vmode,epsilond,T));
        else
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*integral(@(x) mmInt11(q1,q2,lambda,x,epsilond,vmode,T,muL,muR),-Inf,Inf);
        end        
    else % These transfer one electron from lead R to lead L,
        if regularizedON
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*(sumMMr11(q1,q2,lambda,muR,muL,vmode,epsilond,T)+...
                sumMMMMrs11(q1,q2,lambda,muR,muL,vmode,epsilond,T));
        else
            w = s_(n1,n2)/(2*pi*hbar_eV)*gammaL*gammaR*integral(@(x) mmInt11(q1,q2,lambda,x,epsilond,vmode,T,muR,muL),-Inf,Inf);
        end      
    end
end

w(isnan(w)) = 0;

end

function mu = whichMu(muL,muR,lead)
% L(R) with a = +1 (a = -1)
if lead == 1
    mu = muL;
else
    mu = muR;
end

end
function gamma = whichGamma(gammaL,gammaR,lead)
% L(R) with a = +1 (a = -1)
if lead == 1
    gamma = gammaL;
else
    gamma = gammaR;
end

end

function result = mmInt11(q1,q2,lambda,epsilon,epsilond,vmode,T,muL,muR)

maxq = max(q1,q2);
minq = min(q1,q2);
tol = 0;

for q = minq:maxq
    tol = tol + mFCMatrix(q2,q,lambda).*conj(mFCMatrix(q1,q,lambda))./(epsilond-epsilon+(q2-q)*vmode);
end
result = tol.^2.*fermi(epsilon,muL,T).*(1-fermi(epsilon+(q1-q2)*vmode,muR,T));

end

function result = mmInt00(q1,q2,lambda,epsilon,epsilond,vmode,T,muL,muR)

maxq = max(q1,q2);
minq = min(q1,q2);
tol = 0;
for q = minq:maxq
    tol = tol + mFCMatrix(q2,q,lambda).*conj(mFCMatrix(q1,q,lambda))./(epsilon-epsilond+(q1-q)*vmode);
end
result = tol.^2.*fermi(epsilon,muL,T).*(1-fermi(epsilon+(q1-q2)*vmode,muR,T));

end

function s = s_(n1,n2)
% spin degeneracy 
if n1 == 0
    s = 2;
elseif n1 == 1
    if n2 == 0
        s = 1;
    elseif n2 == 1
        s = 2;
    end
end

end
function [db,CommonCommentsValue,flag] = get_set_values(db,name,value)
flag = false; % exist?
if isfield(db,name)
    CommonCommentsValue = find(cellfun(@(x) isequaln(x,value),db.(name)));

    if CommonCommentsValue ~= 0
        flag = true;
    else
        if isempty(CommonCommentsValue)
            CommonCommentsValue = numel(db.(name));
        end
        CommonCommentsValue = CommonCommentsValue+1;
        db.(name){CommonCommentsValue} = value;        
    end
else
    db.(name) = {value};
    CommonCommentsValue = 1;
end
end