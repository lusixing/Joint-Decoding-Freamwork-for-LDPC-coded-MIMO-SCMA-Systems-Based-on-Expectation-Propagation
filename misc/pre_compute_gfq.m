function pre_compute_gfq(q)
 global add_result mul_result div_result

 for x=0:q-1
     for y=0:q-1
         gfadd =gf(x,log2(q))+gf(y,log2(q));
         gfmul =gf(x,log2(q))*gf(y,log2(q));
         add_result(x+1,y+1) =gfadd.x;
         mul_result(x+1,y+1) =gfmul.x;
         
         if y>0
            gfdiv =gf(x,log2(q))/gf(y,log2(q));
            div_result(x+1,y+1) =gfdiv.x;
         end
     end
 end


end