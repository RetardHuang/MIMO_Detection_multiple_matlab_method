%% Reference
% Tan, X., Ueng, Y. L., Zhang, Z., You, X., & Zhang, C.,
% "A Low-Complexity Massive MIMO Detection Based on Approximate Expectation Propagation".
% IEEE Transactions on Vehicular Technology, 68(8),2019. 7260-7272.
function pos_out = EPA(TxRx,N0,Iter,HTH,HTy,H,y,Constellations) 
% inputs: TxRx--MIMO system parameters, N0--noise variance, %Iter--iterations, HTH--real-valued H'*H,
% HTy--real-valued H'*y, H--real-valued channel, y--real-valued received symbols, 
% Constellations--the real-valued constellations(such as [-3 -1 1 3] for 16QAM), 
% output: pos_out--positions for the hard decision constellations 
% Remark: TxRx.Constellations--the complex-valued constellations,
% TxRx.Modulation_order--modulation order(2 for qpsk,4 for 16-QAM,6 for 64-QAM,8 for 256-QAM)

 Es = TxRx.Es/2;  % real value 
 n = 2 * TxRx.Ntx;  %antennas number
 N0 = N0/2;
 lambda = 1/Es;
 belta = 0.8; %the parameter is the damping factor in the reference
 ini = 'init1';
 % initialization
 pos_out = zeros(n/2,Iter);
 mA = 1/N0*(HTH);
 mb = 1/N0*HTy;
 switch(ini)
     case 'init1'
         gama   = zeros(n , 1);
         sigma  = lambda* ones(n , 1);
         Bold_V = ( mA + diag(sigma))\eye(size(mA));  % Eq(7)
         Bold_E = Bold_V*( mb + gama);                % Eq(8)
         Vq     = diag(Bold_V);
         Eq     = Bold_E;
         h      = Vq ./(1 - Vq .* sigma);
         t      = h.*(Eq./Vq - gama);
         V      = diag(mA);
         rou    = t./h;
     case 'init2'
         V      = diag(mA);
         rou    = 0;
         h      = 1./V;
         t      = rou.*h;
     otherwise
         error('detector is not supported.');
 end
 
 for i = 1:Iter
     rou_old = rou;
     value = -0.5*abs(ones(n,1)*Constellations - t*ones(1,sqrt(2^TxRx.Modulation_order))).^2./(h*ones(1,sqrt(2^TxRx.Modulation_order)));
     dex = max(value,[],2);
     prob_muti = exp(value - dex*ones(1,sqrt(2^(TxRx.Modulation_order))));  %calculate the probabilities used h and t
     prob_muti = prob_muti./(sum(prob_muti,2)*ones(1,sqrt(2^TxRx.Modulation_order)));  % normalization
     
     miu = prob_muti * Constellations.';  % Eq(22)
     delta = sum(prob_muti .* abs(ones(2*TxRx.Ntx,1)*Constellations - miu*ones(1,sqrt(2^TxRx.Modulation_order))).^2,2); % Eq(23)
%      delta(delta<1e-4) = 1e-4;
     m = y - H*miu;  % Eq(18a)
     rou = 1/N0*H'*m + V.*miu; % Eq(18b)
     rou  = belta*rou + (1-belta)*rou_old; % Eq(24)
     
     h = 1./V;       % Eq(21)
     t = rou.*h;     % Eq(21)
     
     % turn to complex-domain
     mean = reshape(miu,TxRx.Ntx,[]);
     symb = mean(:,1) + 1j*mean(:,2);
     var = reshape(delta,TxRx.Ntx,[]);
     vv = var(:,1) + var(:,2);
     
     % hard decision(calculate the probabilities of soft decision)
     prob_muti1 = exp(-0.5*abs(ones(TxRx.Ntx,1)*TxRx.Constellations- symb*ones(1,2^(TxRx.Modulation_order))).^2./(vv*ones(1,2^TxRx.Modulation_order)))./(sqrt(2*pi.*(vv*ones(1,2^TxRx.Modulation_order))));
     prob = prob_muti1./(sum(prob_muti1,2)*ones(1,2^TxRx.Modulation_order)); % normalization
     [~,pos_out(:,i)] = max(prob,[],2);
     % early terminated 
%      flag = Vq - Vqf(:,i);
%      if sum(abs(flag))<1e-4
%          pos_out(:,i+1:Iter) = pos_out(:,i)*ones(1,Iter-i);
%          break;         
%      end
 end
end