function [H]=genH_GFq_random(rows,cols,q)

row_flag(1:rows)=0;
parity_check=zeros(rows,cols);

%add bits_per_col 1's to each column with the only constraint being that the 1's should be
%placed in distinct rows
bits_per_col=3;
for i=1:cols
    a=randperm(rows);
    for j=1:bits_per_col
        parity_check(a(j),i)=1;
        row_flag(a(j))=row_flag(a(j))+1;
    end
end

max_ones_per_row=ceil(cols*bits_per_col/rows);

%add 1's to rows having no 1(a redundant row) or only one 1(that bit in the codeword becomes
%zero irrespective of the input)
for i=1:rows
    if row_flag(i)==0    
        for k=1:2
            j=unidrnd(cols);
            while parity_check(i,j)==1
                j=unidrnd(cols);
            end
            parity_check(i,j)=1;        
            row_flag(i)=row_flag(i)+1;  
        end
    end
    if row_flag(i)==1    
        j=unidrnd(cols);
        while parity_check(i,j)==1
            j=unidrnd(cols);
        end
        parity_check(i,j)=1;
        row_flag(i)=row_flag(i)+1;
    end
end

%try to distribute the ones so that the number of ones per row is as uniform as possible
for i=1:rows
    j=1;
    a=randperm(cols);
    while row_flag(i)>max_ones_per_row
        if parity_check(i,a(j))==1 
            newrow=unidrnd(rows);
            k=0;
            while (row_flag(newrow)>=max_ones_per_row | parity_check(newrow,a(j))==1) & k<rows
                newrow=unidrnd(rows);
                k=k+1;
            end
            if parity_check(newrow,a(j))==0
                parity_check(newrow,a(j))=1;
                row_flag(newrow)=row_flag(newrow)+1;
                parity_check(i,a(j))=0;
                row_flag(i)=row_flag(i)-1;
            end
        end
        j=j+1;
    end
end

%try to eliminate cycles of length 4 in the factor graph
for loop=1:10
    chkfinish=1;
    for r=1:rows
        ones_position=find(parity_check(r,:)==1);
        ones_count=length(ones_position);
        for i=[1:r-1 r+1:rows]
            common=0;
            for j=1:ones_count
                if parity_check(i,ones_position(j))==1
                    common=common+1 ;
                    if common==1
                        thecol=ones_position(j);
                    end
                end
                if common==2
                    chkfinish=0; 
                    common=common-1;
                    if (round(rand)==0)         
                        coltoberearranged=thecol;           
                        thecol=ones_position(j);
                    else
                        coltoberearranged=ones_position(j); 
                    end
                    parity_check(i,coltoberearranged)=3; %make this entry 3 so that we dont use
                    %of this entry again while getting rid
                    %of other cylces
                    newrow=unidrnd(rows);
                    iteration=0;     
                    while parity_check(newrow,coltoberearranged)~=0 & iteration<5
                        newrow=unidrnd(rows);
                        iteration=iteration+1;
                    end
                    if iteration>=5  
                        while parity_check(newrow,coltoberearranged)==1
                            newrow=unidrnd(rows);
                        end
                    end
                    parity_check(newrow,coltoberearranged)=1;
                end
            end
        end
    end
    
    if chkfinish
        break
    end
end

%replace the 3's with 0's
parity_check=parity_check==1;
H=zeros(size(parity_check));
%GF(q)
for r=1:rows
    for c=1:cols
        if parity_check(r,c)==1
            H(r,c) =randi([1,q-1]);
        end
    end
end


