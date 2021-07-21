%% 
%  Demo for non-binary ldpc coded sparse code multiple access system with 
%  multiple trceive antennas, working under downlink spatial modulation mode.  
%  Code written by Sixing Lu, Email:lusixingjnu@foxmail.com
%  2021 Mar. 15 
%% 
clear
addpath(genpath('misc'))
loadH = true;   %use the parity check matrix generated last time

Nit = 8;    %number of interation for decoding algorithm
snr_db = 5.5:0.1:6.0;
snr_linear = 10.^(snr_db/10);

n_trans = 1e1;  %number of transmit blocks

Nr =4;  %Number of receive antennas for each user
Nt =4;  %Number of transmit antennas at BS

b_err =zeros(length(snr_db),1);
ber =zeros(length(snr_db),1);

M1 = Nt;
[CB,F,N,V,M2] = get_default_CB();   %N: orthogonal channels V: number of users 

Nb1 = log2(Nt);
Nb2 = log2(M2);

M =M1*M2;

N_msg =256; %message length in M symbols, m_len==rows
b_len =N_msg*log2(M);

rate = 0.5; %chennel coding rate
rows = N_msg/rate - N_msg;
cols = N_msg + rows;
q = M;
N_scma = cols;

global add_result mul_result div_result
add_result = zeros(q,q);
mul_result = zeros(q,q);
div_result = zeros(q,q);
pre_compute_gfq(q);

if ~(loadH && exist("sim_data.mat","file"))
    H = genH_GFq_random(rows,cols,M);   %generate parity check matrix over gf(q)
    [P,shifted_cols] =H2P_GFq(H,q);    %get generator matrix from H
else
    load("sim_data","H","P","shifted_cols")
end

Nc = N*Nr;

for sj=1:length(snr_linear)
    snr = snr_linear(sj);
    for nf=1:n_trans
        b1 = randi([0,1],[Nb1,N_msg,V]); %data stream 1 in bits, mapping the antennas index
        b2 = randi([0,1],[Nb2,N_msg,V]); %data stream 2 in bits, mapping the codework   
        
        s1 = zeros(N_msg,V);
        s2 = zeros(N_msg,V);
        s = zeros(N_msg,V);
        
        h = randn(N,Nr,Nt) + 1i*randn(N,Nr,Nt); %channel
        
        for i=1:N_msg
            for v=1:V
                s1(i,v) = bi2de(b1(:,i,v)','left-msb'); %antennas index
                s2(i,v) = bi2de(b2(:,i,v)','left-msb');
                s(i,v) = (s1(i,v))*M2 + s2(i,v);
            end
        end
        %% construct the unified codebook and factor graph 
        CB_V = zeros(Nc,M,V);
        for v=1:V
            for m1=1:M1
                for m2=1:M2
                    m = (m1-1)*M2 +m2;
                    for nr=1:Nr
                        CB_V((nr-1)*N+1:nr*N,m,v) = h(:,nr,m1).*CB(:,m2,v);
                    end
                end
            end
        end
        
        Fc = zeros(Nc,V);
        
        for nr=1:Nr
            Fc((nr-1)*N+1:nr*N,:) =F;
        end
        %% ----------------NB-ldpc encoding---------
        y = zeros(Nc,N_scma);
        C1 =zeros(cols,V);
        C2 =zeros(cols,V);
        C =zeros(cols,V);
        X = zeros(N,Nt,N_scma);
        
        for v=1:V
            C(:,v) =ldpc_encode_GFq(s(:,v)',P,q,shifted_cols);   %ldpc encoding over gf(q)
            %Parity_equation_check(q,H,C(:,v))
        end
       %%  
        for v=1:V
            for c=1:cols
                C1(c,v) = floor(C(c,v)/M2);
                C2(c,v) = C(c,v) - M2*C1(c,v);
                X(:,C1(c,v)+1,c) = X(:,C1(c,v)+1,c) + CB(:,C2(c,v)+1,v); %transmit codeword using spatial modulation mode
            end
        end
        
        for c=1:N_scma
            for nr=1:Nr
                for nt=1:Nt
                    y(N*(nr-1)+1:N*nr,c) = y(N*(nr-1)+1:N*nr,c) + h(:,nr,nt).*X(:,nt,c);
                end
            end

            noise_power = norm(y(:,c),'fro')/sqrt(snr);
            noise = randn(Nc,1) + 1i*randn(Nc,1);
            noise = noise_power*noise./norm(noise,'fro');    
            
            y(:,c) = y(:,c) + noise;
        end
        
        %--------------decoding--------------------------------------
        [s_dec,b_dec] = joint_scma_nb_ldpc_ep_decode_VCB_test1(y, CB_V, Fc, H, shifted_cols, snr, Nit);
        %------------------------------------------------------------
        for v=1:V
            b_err(sj) =b_err(sj) + sum(sum(b_dec(:,:,v)~= [b1(:,:,v);b2(:,:,v)]));
        end
        fprintf("simulation progress: %.2f %%\n",100*((sj-1)*n_trans + nf )/(n_trans*length(snr_db)))
    end
end


for i=1:length(snr_db)
       ber(i) = b_err(i)/(V*(Nb1+Nb2)*N_msg*n_trans);
end


ber_treshold = 1/(V*(Nb1+Nb2)*N_msg*n_trans);
ber = max(ber,ber_treshold);

save("sim_data")


