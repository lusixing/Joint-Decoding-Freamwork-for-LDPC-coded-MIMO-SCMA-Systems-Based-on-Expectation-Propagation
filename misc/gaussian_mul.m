function [m_new,v_new] = gaussian_mul(m1,v1,m2,v2)
    
   v_new = 1/(1/v1 +1/v2);
   m_new = v_new *(m1/v1 +m2/v2);

end
