function dF = BW(y,F,Nnodes)

beta = 1/0.65; 
gamma = 1/0.41; 
tau = 0.98; 
alpha = 0.32;
p = 0.34;

dF = zeros(Nnodes,4);  

dF(:,1) = y - beta*F(:,1) - gamma*(F(:,2)-1); 
dF(:,2) = F(:,1);                             
dF(:,3) = 1/tau*(F(:,2)-F(:,3).^(1/alpha));   
dF(:,4) = 1/tau*(F(:,2)/p.*(1-(1-p).^(1./F(:,2)))-F(:,4)./F(:,3).*F(:,3).^(1/alpha)); %dq: deoxyHb

end


