function C = matrix_mul_gfq(q,A,B)
global add_result mul_result

[a_rows,a_cols] = size(A);
[b_rows,b_cols] = size(B);

C =zeros(a_rows,b_cols);
for i =1:a_rows
    for j =1:b_cols
        temp =0;
        for l=1:a_cols
            if A(i,l)>=q || B(l,j)>=q
                ME = MException('GfValErr:illegalval', ...
                    'illegal value for gf(q)',q);
                throw(ME)
            end
            temp =add_result(temp+1,mul_result(A(i,l)+1 ,B(l,j)+1)+1 );
            
        end
        C(i,j)=temp;
    end
end

