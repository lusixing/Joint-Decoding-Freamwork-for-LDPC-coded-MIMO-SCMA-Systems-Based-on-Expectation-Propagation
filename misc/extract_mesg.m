function [u]= extract_mesg(c,shifted_cols)
rows=length(shifted_cols);
for i=1:rows
   if shifted_cols(i)~=0
      temp=c(i);
      c(i)=c(shifted_cols(i));
      c(shifted_cols(i))=temp;
   end
end
cols=length(c);
u=c(rows+1:cols);
