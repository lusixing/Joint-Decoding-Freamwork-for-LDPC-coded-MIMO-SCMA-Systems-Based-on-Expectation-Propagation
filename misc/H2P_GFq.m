%compute the generator matrix from parity check matrix over gf(q)
function [P,rearranged_cols]=H2P_GFq(H,q)
rows =size(H,1);
cols =size(H,2);

rearranged_cols =zeros(1, rows);

%% gaussian elimination over gf(q)
for r =1:rows
    if H(r,r)==0 %|| sum(H(1:r-1,r))~=0
       for c=r+1:cols
           if H(r,c)~=0 %&& sum(H(1:r-1,c))==0
               rearranged_cols(r) =c;
               temp =H(:,r);
               H(:,r)=H(:,c);
               H(:,c)=temp;
               break;
           end
           if c==cols 
                ME = MException('HErr:illegalH', ...
                    'illegal parity check matrix H',H);
                throw(ME)
           end
       end 
    end
    H(r,:) =div_gfq(q,H(r,:),H(r,r))';
    
    for r2=1:rows
       if r2~=r && H(r2,r)~=0
           %H(r2,:) =H(r2,:) - H(r2,r).*H(r,:);
           H(r2,:) = add_gfq(q,H(r2,:),mul_gfq(q, H(r2,r),H(r,:))')';
       end
    end

end

P =H(:,rows+1:cols);
