clc;
clear all;

load em_SC.mat
load em_FC.mat

FC(FC < 0) = 0; %Convert the negative value of the FC matrix to 0

%Normalize the SC
SC = SC./max(max(SC));
SC = log(SC+1);

%find out number of brain regions
NumC = length(diag(SC)); 

FC_mask = tril(ones(size(FC,1),size(FC,1)),0);
y = FC(~FC_mask); 
T = length(y);    

len=29;
G=linspace(0.1,1.5,len);
state_corr=zeros(1,len);

 for i=1:len
[FC_sim, CC_check, BOLD, y_neuro_cut_E, y_neuro_cut_I, H_neuro_cut_E, H_neuro_cut_I, FC_simR] = DMF_E_I_eul(G(i),SC,y,FC_mask,18,0.72);
   state_corr(i)=CC_check;
 end

 %fig2b
 figure
 plot(G,state_corr);


