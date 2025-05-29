function [dsn, dsg, rn, rg]=E_I_DMF(G,sn,sg,SC)

taon=0.1; 
taog=0.01; 
gamma=0.641;
J_E=0.15; 
J_I=1; 
I0=0.382; 
Jexte=1; 
Jexti=0.7; 
W_EE=0.1; 
W_EI=1;
W_IE=1; 
W_II=1; 

ae=310; 
be=125;
de=0.16;

ai=615;
bi=177;
di=0.087;

xn=I0*Jexte+W_EE*J_E*sn+G*J_E*SC*sn-W_IE*J_I*sg;
xg=I0*Jexti+W_EI*J_E*sn-W_II*J_I*sg;

He=(ae*xn-be)./(1-exp(-de*(ae*xn-be)));
Hi=(ai*xg-bi)./(1-exp(-di*(ai*xg-bi)));

rn=He;
rg=Hi;
dsn=-sn/taon+(1-sn)*gamma.*rn;    
dsg=-sg/taog+rg;

end



