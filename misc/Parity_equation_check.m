function bool =Parity_equation_check(q,H,u)
rows =size(H,1);
check =gf(zeros(1,rows),log2(q));
H = gf(H,log2(q));
u = gf(u,log2(q));

for r=1:rows
   ind=find(H(r,:)~=0);
   check(r) =gf(0,log2(q));
   for i=1:length(ind)
        check(r) =check(r) +u(ind(i))*H(r,ind(i)); 
   end
   if check(r)~=gf(0,log2(q))
       bool =0;
       return
   end
end

bool=1;