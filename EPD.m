%% Reference
% J. Céspedes, P. M. Olmos, M. Sanchez-Fernandez, and F. Perez-Cruz,
% “Expectation propagation detection for high-order high-dimensional
% MIMO systems,” IEEE Trans. Commun., vol. 62, no. 8, pp. 2840–2849,
% Aug. 2014.

function pos_out = EPD(TxRx,N0,Iter,HTH,HTy,Constellations,symbol) 
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
sigma = 1/Es* ones(n , 1);
% sigma = 0.001* ones(n , 1);

efcelong = 5e-7;
belta = 0.2;


 mA = 1/N0*(HTH);
 mb = 1/N0*HTy;

 Bold_V = ( mA + diag(sigma))\eye(size(mA));  % Eq(28)
 Bold_E = Bold_V*( mb + gama);                % Eq(29)
 Vq = diag(Bold_V);
 Eq = Bold_E;
%   Vq = sum(Vq)/length(Vq)*ones(length(Vq),1);
 for i = 1:Iter 
     
     Eqf(:,i) = Eq;
     Vqf(:,i) = Vq;
     sigma_re = sigma;
     gama_re =gama;

     %% main updating

     h = Vq ./(1 - Vq .* sigma);  % Eq(31)
     t = h.*(Eq./Vq - gama);      % Eq(32)
     
     % calculate the sampling probabilities of cavity distributions
     value = -0.5*abs(ones(n,1)*Constellations - t*ones(1,sqrt(2^TxRx.Modulation_order))).^2./(h*ones(1,sqrt(2^TxRx.Modulation_order)));
     dex = max(value,[],2);
     prob_muti = exp(value-dex*ones(1,sqrt(2^(TxRx.Modulation_order))));
     prob_muti = prob_muti./(sum(prob_muti,2)*ones(1,sqrt(2^TxRx.Modulation_order)));
     
     % step 2) and Eq(33)
     miu  = prob_muti * Constellations.';
     delta = sum(prob_muti .* abs(ones(n,1)*Constellations - miu*ones(1,sqrt(2^TxRx.Modulation_order))).^2,2);
     delta = max(efcelong ,delta);
%      delta = sum(delta)/length(delta)*ones(length(delta),1);
     % when delta(i)>h(i) is true, the present updating is canceled. 
     gama = gama_re + belta*(miu./delta - Eq./Vq).*(delta<h);  % Eq(37),Eq(31)and Eq(32) are uesed here.
     sigma = sigma_re + belta*(1./delta - 1./Vq).*(delta<h);   % Eq(38)
%      gama = miu./delta - t./h;  % Eq(37),Eq(31)and Eq(32) are uesed here.
%      sigma = 1./delta - 1./h;   % Eq(38)
%      sigma = sum(sigma)/length(sigma)*ones(length(sigma),1);
     % Eq(28) and Eq(29) for next iteration
     Bold_V = ( mA + diag(sigma))\eye(size(mA)); %main complexity load (iterative matrix inversion)
     Bold_E = Bold_V*( mb + gama);
     Vq = diag(Bold_V);
%      Vq = sum(Vq)/length(Vq)*ones(length(Vq),1);
     Eq = Bold_E;
    %%
%      [f,xi] = ksdensity(Vqf(:,i) - Vq);
%      
%      plot(xi,f);
%      hold on
     % turn to complex-domain
     mean = reshape(Eq,n/2,[]);
     symb = mean(:,1) + 1j*mean(:,2);
%      var = reshape(Vq,n/2,[]);
%      vv = var(:,1) + var(:,2);
%      
%      % hard decision(calculate the probabilities of soft decision)
%      value = -0.5*abs(ones(n/2,1)*TxRx.Constellations- symb*ones(1,2^(TxRx.Modulation_order))).^2./(vv*ones(1,2^TxRx.Modulation_order));
%      dex = max(value,[],2);
%      prob_muti1 = exp(value-dex*ones(1,2^(TxRx.Modulation_order)))./(sqrt(2*pi.*(vv*ones(1,2^TxRx.Modulation_order))));
%      prob = prob_muti1./(sum(prob_muti1,2)*ones(1,2^TxRx.Modulation_order));
%      % decision for every iteration, if only the final results is needed, the calculation can finished outside the loop
%      [~,pos_out(:,i)] = max(prob,[],2); 
     [~,pos_out(:,i)] = min(abs(symb - TxRx.Constellations),[],2);
     % early terminated 
     flag = Vq - Vqf(:,i);
     if sum(abs(flag))<1e-4
         pos_out(:,i+1:Iter) = pos_out(:,i)*ones(1,Iter-i);
         break;         
     end

 end %iteration end
% close all
% plot(sum(Vqf)/n)

end %function end

