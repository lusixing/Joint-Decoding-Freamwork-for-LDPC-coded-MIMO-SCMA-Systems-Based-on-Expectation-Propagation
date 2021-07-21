%compute the log of complex gaussian distribution pdf
function p = cnormpdf_l(z,m,v)
p = -conj(z-m)*(z-m)/v;