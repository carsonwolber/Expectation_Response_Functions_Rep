function out = gomega(w,A,B,v,I)

T1 = A/(I-B.*exp(-1i.*w));
T2 = A/(I-B.*exp(1i.*w));
out = diag(1/(2*pi)*(T1*v*(T2*v).')./((1-exp(-1i.*w)).*(1-exp(1i.*w))));