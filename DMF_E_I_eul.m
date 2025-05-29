function [FC_sim, CC_check, BOLD_d, y_neuro_cut_E, y_neuro_cut_I, H_neuro_cut_E, H_neuro_cut_I, FC_simR] = DMF_E_I_eul(G,SC,y_FC,FC_mask,Tepochlong,TBOLD)

kstart = 0;  
Tpre = 60*2; 
kend = Tpre+60*Tepochlong; 

dt_l = 0.005;    

dt = dt_l;  
dtt = 0.005; 

k_P = kstart:dt:kend;   
k_PP = kstart:dtt:kend;

Nnodes = size(SC,1);
Bsamples = length(k_PP);

ysn = zeros(Nnodes,1);
ysg = zeros(Nnodes,1);

zT = zeros(Nnodes,Bsamples);
fT = zeros(Nnodes,Bsamples);
fT(:,1) = 1;
vT = zeros(Nnodes,Bsamples);
vT(:,1) = 1;
qT = zeros(Nnodes,Bsamples);
qT(:,1) = 1;

F = [zT(:,1) fT(:,1) vT(:,1) qT(:,1)];
ysn(:,1) = 0.001;
ysg(:,1) = 0.001;

sigma = 0.01; 
w_L = length(k_P);
dW = sqrt(dt)*randn(Nnodes,w_L); 
j = 0;

for i = 1:length(k_P)
        
        [dsn, dsg, rn, rg]=E_I_DMF(G,ysn,ysg,SC);
  ysn=ysn+dt*dsn+sigma*dW(:,i);           
  ysg=ysg+dt*dsg+sigma*dW(:,i);

        if mod(i,dtt/dt) == 0
            j = j+1;
            y_neuro_E(:,j) = ysn; 
            H_neuro_E(:,j) = rn;  

            y_neuro_I(:,j) = ysg; 
            H_neuro_I(:,j) = rg;  
        end
        
end

for i = 2:length(k_PP)

            dF = BW(y_neuro_E(:,i-1),F,Nnodes);
            F = F + dF*dtt;
            zT(:,i) = F(:,1);
            fT(:,i) = F(:,2);
            vT(:,i) = F(:,3);
            qT(:,i) = F(:,4);

end

p = 0.34; 
v0 = 0.02;


k1 = 4.1;  
k2 = 0.58;  
k3 =0.53;   


y_BOLD = 100/p*v0*( k1*(1-qT) + k2*(1-qT./vT) + k3*(1-vT) );    

Time = k_PP;

cut_indx = find(Time == Tpre);
BOLD_cut = y_BOLD(:,cut_indx:end);
y_neuro_cut_E = y_neuro_E(:,cut_indx:end);
y_neuro_cut_I = y_neuro_I(:,cut_indx:end);
H_neuro_cut_E = H_neuro_E(:,cut_indx:end);
H_neuro_cut_I = H_neuro_I(:,cut_indx:end);

BOLD_d = simBOLD_downsampling(BOLD_cut,TBOLD/dtt); %down sample 

FC_sim = corr(BOLD_d');
FC_sim=atanh(FC_sim);
FC_sim(FC_sim < 0) = 0; 
FC_simR = FC_sim(~FC_mask);

CC_check = corr(FC_simR,y_FC);

end

