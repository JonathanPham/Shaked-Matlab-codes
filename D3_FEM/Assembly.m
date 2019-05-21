function [K F]=Assembly() %done
global nEquations nElements nDoF ID LM h;
K = spalloc(nEquations,nEquations,20*nEquations);
F = spalloc(nEquations,1,nEquations);
for e=1:nElements
    [k_e, f_e, f_g, f_h]=linhex(e); %just need to define this
    P=LM(:,e);
    I=(P~=0);
    P=P(I);
    K(P,P)= K(P,P)+k_e(I,I);
    F(P)=F(P)+f_e(I) + f_g(I) +f_h(I);
end
end