function [u]=ldpc_encode_GFq(s,P,q,rearranged_cols)

s = reshape(s,[length(s),1]);
c =matrix_mul_gfq(q,P,s);

u=[c' s']; 

u =reorder_pos(u,rearranged_cols);
