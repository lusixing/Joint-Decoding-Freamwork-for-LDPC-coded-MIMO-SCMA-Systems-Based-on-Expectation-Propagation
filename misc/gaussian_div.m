function [m_new,v_new] = gaussian_div(m1,v1,m2,v2)

v_new = 1/(1/v1 -1/v2);

if isnan(v_new) || ~isfinite(v_new)
    ME = MException('NumericalErr:badval', ...
        'value of %s corrupted, check the configuration',v_new);
    throw(ME)
end

if v_new<0
    v_new=1e10;
end

if abs(v_new)<1e-20
    v_new = 1e-20;
end


m_new = v_new*(m1/v1 -m2/v2);
end