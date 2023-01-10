%% Reference
% J. Céspedes, P. M. Olmos, M. Sanchez-Fernandez, and F. Perez-Cruz,
% “Expectation propagation detection for high-order high-dimensional
% MIMO systems,” IEEE Trans. Commun., vol. 62, no. 8, pp. 2840–2849,
% Aug. 2014.

function pos_out = EPL(TxRx,N0,Iter,H,y,Constellations,symbol) 
%EPD 此处显示有关此函数的摘要
% inputs: TxRx--MIMO system parameters, N0--noise variance, %Iter--iterations, HTH--real-valued H'*H,
% HTy--real-valued H'*y, Constellations--the real-valued constellations(such as [-3 -1 1 3] for 16QAM),
% symbols--for debug
% output: pos_out--positions for the hard decision constellations 

sy=TxRx.Constellations(symbol).';a = [real(sy);imag(sy)]; %just for debug

Es = sum(abs(TxRx.Constellations).^2)/(2^TxRx.Modulation_order); % average energy for transmitted symbols
% Es = Es/2;  % real value 
n = 2 * TxRx.Ntx;  %antennas number
% N0 = N0/2;

% initialization
pos_out = zeros(n/2,Iter);

gama = zeros(n , 1);
% sigma = 1/Es* ones(n , 1);
sigma = 1/Es;

beta = 0.2;
efcelong = 5e-5;

[~,dia,V] = svd(H);

% HTH = H'*H;
% [dia1,V1] = Jaco(HTH);
HTy = H'*y;
mb = 1/N0*HTy;

 Inv = N0./((diag(dia)).^2 + sigma*N0);
%  Inv1 = N0./(diag(dia1) + sigma*N0);
 Ve = ( mb + gama);
 Ve1 = V'*Ve;
 Ve2 = Inv.*Ve1;
 Eq = V*Ve2;  % Eq(28)
%  Vq = sum(Inv)/length(Inv)*ones(length(Eq),1);
 Vq = sum(Inv)/length(Inv);
%  Vf1 =  V1'*Ve;
%  Vf2 = Inv1.*Vf1;
%  Eq1 = V1*Vf2; 
 for i = 1:Iter 
     
     Vqf(i) = Vq;
     sigma_re = sigma;
     gama_re =gama;

     %% main updating
% 
%      h = Vq ./(1 - Vq .* sigma);  % Eq(31)
%      t = h.*(Eq./Vq - gama);      % Eq(32)
     h = Vq /(1 - Vq * sigma);  % Eq(31)
     t = h*(Eq/Vq - gama);      % Eq(32)
     
     % calculate the sampling probabilities of cavity distributions
%      value = -0.5*abs(ones(n,1)*Constellations - t*ones(1,sqrt(2^TxRx.Modulation_order))).^2./(h*ones(1,sqrt(2^TxRx.Modulation_order)));
     value = -0.5*abs(ones(n,1)*Constellations - t*ones(1,sqrt(2^TxRx.Modulation_order))).^2./h;
     dex = max(value,[],2);
     prob_muti = exp(value-dex*ones(1,sqrt(2^(TxRx.Modulation_order))));
     prob_muti = prob_muti./(sum(prob_muti,2)*ones(1,sqrt(2^TxRx.Modulation_order)));
     
     % step 2) and Eq(33)
     miu  = prob_muti * Constellations.';
     delta = sum(prob_muti .* abs(ones(n,1)*Constellations - miu*ones(1,sqrt(2^TxRx.Modulation_order))).^2,2);
     delta1 = max(efcelong ,delta);
%       delta = sum(delta1)/length(delta1)*ones(length(delta1),1);
     delta = sum(delta1)/length(delta1);
     % when delta(i)>h(i) is true, the present updating is canceled. 
%      gama = gama_re + beta*(miu/delta - Eq/Vq);  % Eq(37),Eq(31)and Eq(32) are uesed here.
%      sigma = sigma_re + beta*(1/delta - 1/Vq);   % Eq(38)
     gama = miu/delta - t/h;  % Eq(37),Eq(31)and Eq(32) are uesed here.
     sigma = 1/delta - 1/h;   % Eq(38)
     
     % Eq(28) and Eq(29) for next iteration
      Inv = N0./((diag(dia)).^2 + sigma*N0);
      Ve = ( mb + gama);
      Ve1 = V'*Ve;
      Ve2 = Inv.*Ve1;
      Eq = V*Ve2;  % Eq(28)
      Vq = sum(Inv)/length(Inv);

     % turn to complex-domain
     mean = reshape(Eq,n/2,[]);
     symb = mean(:,1) + 1j*mean(:,2);

     [~,pos_out(:,i)] = min(abs(symb - TxRx.Constellations),[],2);
     % early terminated 
     flag = Vq - Vqf(i);
     if sum(abs(flag))<1e-4
         pos_out(:,i+1:Iter) = pos_out(:,i)*ones(1,Iter-i);
         break;         
     end

 end %iteration end
% close all
% plot(sum(Vqf)/n)

end %function end

