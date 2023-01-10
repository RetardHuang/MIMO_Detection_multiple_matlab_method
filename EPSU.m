function pos_out = EPSU( TxRx,N0,Iter,HTH,HTy,Constellations)
%EPNEW 此处显示有关此函数的摘要
%   此处显示详细说明
pos_out = zeros(TxRx.Ntx,Iter);
Es = sum(abs(TxRx.Constellations).^2)/(2^TxRx.Modulation_order);
NN = 2 * TxRx.Ntx;
% load('sym_pos.mat');
% Constellations = [-3 -1 1 3];
% proflag = ones(NN,length(Constellations));
 
 gama = zeros(NN , 1);
 sigma = 1/Es * ones(NN , 1);

 
 h = zeros(NN , 1);
 t = zeros(NN , 1);
 
 miu = zeros(NN , 1);
 delta = zeros(NN , 1);
 
 prob_muti = zeros(NN ,TxRx.Modulation_order);

 efcelong = 5e-7;
 belta = 0.2;
 alpha = 0.5;
 
 mA = 1/N0*(HTH);
 mb = 1/N0*HTy;
 
 Bold_V = ( mA + diag(sigma))\eye(size(mA));
 Bold_E = Bold_V*( mb + gama);
 Vq = real(diag(Bold_V));
 Eq = Bold_E;
 
 sigma_old = sigma;
 for i = 1:Iter        
     Eqf(:,i) = Eq;
     Vqf(:,i) = Vq;
     for n = 1:NN
%          
         h(n) = Vq(n)/(1 - Vq(n) * sigma(n));
         t(n) = h(n)*(Eq(n)/Vq(n) - gama(n));
         prob_muti(n,:) = exp(-alpha*abs(Constellations - t(n)).^2/h(n));
         prob_muti(n,:) = prob_muti(n,:)/sum(prob_muti(n,:));

         miu(n)  = sum(prob_muti(n,:) .* Constellations);
         delta(n) = sum(prob_muti(n,:) .* abs(Constellations - miu(n)).^2);
         delta(n) = max(efcelong ,delta(n));

         sigma(n) = belta*(1/delta(n) - 1/h(n)) + (1 - belta)*sigma(n);
         gama(n) = belta*(miu(n)/delta(n) - t(n)/h(n)) + (1 - belta)*gama(n);
        
%          Bold_V =  Bold_V - ( Bold_V * diag(sigma - sigma_old)*Bold_V)/(1+ Bold_V(n,n)*sum(sigma - sigma_old));
%          deltalamda = sigma(n) - sigma_old(n);
%          para = deltalamda/(1+deltalamda*Vq(n));
%          Bold_V =  Bold_V - para*(Bold_V(:,n)* Bold_V(n,:));
         Bold_V = ( mA + diag(sigma))\eye(size(mA));
         Bold_E = Bold_V*(mb + gama);
         Vq = real(diag(Bold_V));
         Eq = Bold_E;
         sigma_old = sigma;
         
     end
     
     mean = reshape(Eq,TxRx.Ntx,[]);
     symb = mean(:,1) + 1j*mean(:,2);

     [~,pos_out(:,i)] = min(abs(symb - TxRx.Constellations),[],2);
     % early terminated 
     flag = Vq - Vqf(:,i);
     if sum(abs(flag))<1e-4
         pos_out(:,i+1:Iter) = pos_out(:,i)*ones(1,Iter-i);
         break;         
     end
 end

end

