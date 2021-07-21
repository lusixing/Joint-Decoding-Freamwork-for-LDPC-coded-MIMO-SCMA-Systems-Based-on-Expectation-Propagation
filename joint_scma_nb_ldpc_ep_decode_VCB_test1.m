function [s_dec,b_dec] = joint_scma_nb_ldpc_ep_decode_VCB_test1(y,CB_V,Fc,H,shifted_cols,snr,Nit)
%%
% Joint decoding algorithm for NB-LDPC coded SCMA based on expectation
% propagation, using unified codebook and unified factor graph  repersentation.
% Code written by Sixing Lu, Email: lusixingjnu@foxmail.com
%%
% y: received signal
% CB_V: unified codebook
% H: NB-LDPC check matrix over gf(q)
% Fc: unified factor graph
% snr: linear snr
% Nit: Number of iterations
%%

[Nc,N_scma] =size(y);
[~,M,V] = size(CB_V);
[rows,cols] = size(H);

q =M;
N_msg= cols-rows;

s_dec = zeros(N_msg,V);
b_dec = zeros(log2(M),N_msg,V);

mu_n_2_v = zeros(Nc,V,N_scma);
var_n_2_v =100*ones(Nc,V,N_scma);
mu_v_2_n = zeros(Nc,V,N_scma);
var_v_2_n =100*ones(Nc,V,N_scma);

mu_vn = zeros(Nc,V,N_scma);
var_vn = zeros(Nc,V,N_scma);

q_xy = ones(V,M,N_scma);
q_xy_l = zeros(V,M,N_scma);

Ifc = 1/q*ones([q,size(H),V]);  %gf(q)
Icf = zeros([q,size(H),V]);

c_dec =zeros(N_scma,V);

Ind_c = cell(cols,1);
Ind_r = cell(rows,1);

Ind_v = cell(V,1);
Ind_n = cell(Nc,1);

for c=1:cols
    Ind_c{c} = find(H(:,c)~=0);
end

for r=1:rows
    Ind_r{r} = find(H(r,:)~=0);
end

for v=1:V
    Ind_v{v} = find(Fc(:,v)==1);
end

for n=1:Nc
    Ind_n{n} = find(Fc(n,:)==1);
end

for iter =1:Nit
    %% step 1
    for c=1:cols
        noise_pw = norm(y(:,c),'fro')/sqrt(snr);
        var_noise = noise_pw/Nc;
        for n=1:Nc
            idx_n = Ind_n{n};
            for i=1:length(idx_n)
                m_temp =0;
                v_temp =var_noise;
                for j=1:length(idx_n)
                    if j~=i
                        m_temp =m_temp +  mu_v_2_n(n,idx_n(j),c);
                        v_temp =v_temp +  var_v_2_n(n,idx_n(j),c);
                    end
                end
                mu_n_2_v(n,idx_n(i),c) = (y(n,c)-m_temp);
                var_n_2_v(n,idx_n(i),c) =v_temp;
            end
            
        end
        
    end
    %% step 2
    for c=1:cols
        for v=1:V
            idx_v=Ind_v{v};
            qv_temp_l = zeros(1,M);
            for i=1:length(idx_v)
                for m=1:M
                    qv_temp_l(m) = qv_temp_l(m) + cnormpdf_l(CB_V(idx_v(i),m,v),mu_n_2_v(idx_v(i),v,c),var_n_2_v(idx_v(i),v,c));
                end
            end
            
            q_xy_l(v,:,c) = qv_temp_l -log_sum_exp(qv_temp_l);
            q_xy_l(v,:,c) = max(q_xy_l(v,:,c),-100);
            q_xy(v,:,c) =exp(q_xy_l(v,:,c));
        end
    end
    %% step 3 NB-LDPC steps
    % Icf update---------------
    for v=1:V
        for c=1:cols
            ind_c=Ind_c{c};
            for i=1:length(ind_c)
                p_temp = q_xy(v,:,c)'.*prod(Ifc(:,ind_c(ind_c~=ind_c(i)),c,v),2);
                p_temp = p_temp/sum(p_temp);       %normalize
                Icf(:,ind_c(i),c,v) =p_temp;
            end
        end
    end
    
    % Ifc update------------------
    for v=1:V
        for r=1:rows
            ind = Ind_r{r};
            dr =length(ind);
            perm =zeros(dr,q);
            
            for i=1:dr
                perm (i,:) =div_gfq(q,[0:q-1],H(r,ind(i)))+1;
            end
            
            fh = squeeze(Icf(:,r,ind,v))';              %size(Icf)=([q,size(H),J])
            W = hadamard(q);
            Fh =zeros(size(fh));
            
            for i=1:dr
                Fh(i,:) =fh(i, perm(i,:))*W;
            end
            
            for i=1:dr
                temp =zeros(1,q);
                for q_temp=1:q
                    temp(q_temp) =prod(Fh([[1:i-1] [i+1:dr]],q_temp));
                end
                Py_dr =1/q*temp*W;
                Ifc(perm(i,:),r,ind(i),v) =Py_dr';
                
                for m=1:M
                    Ifc(m,r,ind(i),v) = max(Ifc(m,r,ind(i),v),1e-20);
                end
            end
        end
    end
    
    %% step 4
    for c=1:cols
        ind_c =Ind_c{c};
        for v=1:V
            idx_v = find(Fc(:,v)==1);
            mu_vn_temp = zeros(Nc,1);  %complex 
            var_vn_temp = zeros(Nc,1);   %real 
            
            q_c =  prod(Ifc(:,ind_c,c,v),2);
            q_vc = q_xy(v,:,c).*q_c';

            q_vc = q_vc./sum(q_vc);
            
            for i=1:length(idx_v)
                for m=1:M
                    mu_vn_temp(idx_v(i)) = mu_vn_temp(idx_v(i)) + CB_V(idx_v(i),m,v)*q_vc(m);
                end
                for m=1:M
                    var_vn_temp(idx_v(i)) = var_vn_temp(idx_v(i))+ abs(CB_V(idx_v(i),m,v) - mu_vn_temp(idx_v(i)))^2 *q_vc(m);
                    var_vn_temp(idx_v(i)) =max(var_vn_temp(idx_v(i)),1e-10);
                end
            end
            mu_vn(:,v,c) =mu_vn_temp;
            var_vn(:,v,c) =var_vn_temp;
            
            if isnan(sum(mu_vn(:,v,c)))
                ME = MException('NumericalErr:badval', ...
                    'value of %s corrupted, check the configuration',mu_vn);
                throw(ME)
            end
        end
    end
    %toc
    %% step 5
    for c=1:cols
        for v=1:V
            idx_v = Ind_v{v};
            for i=1:length(idx_v)
                [m_temp,v_temp] = gaussian_div(mu_vn(idx_v(i),v,c),var_vn(idx_v(i),v,c),mu_n_2_v(idx_v(i),v,c),var_n_2_v(idx_v(i),v,c));
                mu_v_2_n(idx_v(i),v,c) =m_temp;
                var_v_2_n(idx_v(i),v,c) =v_temp;
            end
        end
    end
end

%% final decision
for c=1:cols
    ind_c=Ind_c{c};
    for v=1:V
        q_c =  prod(Ifc(:,ind_c,c,v),2);
        q_vc = q_xy(v,:,c).*q_c';
        [~,maxidx]=max(q_vc);
        c_dec(c,v)=maxidx-1;
    end
end


for v=1:V
    s_dec(:,v) =extract_mesg(c_dec(:,v),shifted_cols);
    for n=1:N_msg
         b_dec(:,n,v) = de2bi(s_dec(n,v),log2(M),'left-msb');
    end
end

