function J = regularizedJ11(E1,E2,epsilon,T)
beta = 1/(KBoltzmann_ev*T);
[nxep,nyep] = size(epsilon);
ne2 = numel(E2);
if nxep > 1
    J = zeros([ne2,nyep]);
    for jj = 1:nyep
        temp = epsilon(:,jj);
        J(:,jj) = single_regularizedJ(temp');
    end
else
    J = single_regularizedJ(epsilon);
end
    function single_J = single_regularizedJ(epsilon0)
        a1 = 1/2 + 1i*beta*(E2-epsilon0)/(2*pi);
        a2 = 1/2 + 1i*beta*(E1-epsilon0)/(2*pi);
        t1 = digammaFcn(1,a1);
        t2 = digammaFcn(1,a2);
        single_J = beta/(2*pi)*BoseFcn(E2-E1,T).*imag(t1-t2);
        if ~iscolumn(single_J)
            single_J = single_J';
        end
    end
if iscolumn(J)
    J = J';
end
end