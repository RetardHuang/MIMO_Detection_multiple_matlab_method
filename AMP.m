
function pos_out = AMP(TxRx,y,H,N0,Iter,Constellations,pos_in)
%   AMP detector in Massive MIMO
 sym = TxRx.Constellations(pos_in);
 x = [real(sym.');imag(sym.')];%分开
 m = 2*TxRx.Nrx;%接收天线数
 n = 2*TxRx.Ntx;%发射天线数
 beta = n/m;
 sigmas2 = sum(Constellations.^2)/length(Constellations);
    r = y;
    s = zeros(n,1);
    tau_old = beta/N0*sigmas2;   %initial estimation variance
    for t=1:Iter
        z = s + H'*r;
        value = -1/2*(Constellations-z).^2/(N0*(1+tau_old));
        dex = max(value,[],2);
        prob = exp(value - dex);
        prob = prob./sum(prob,2);
        if(tau_old<1e-5)
            break;
        end
        s = prob*Constellations';%星座点的概率加权
        tau_new = beta/N0*(sum(prob*(Constellations').^2 - s.^2)/n);
        r = y - H*s + tau_new/(1+tau_old)*r;
        tau_old = tau_new;
    end
     mean = reshape(s,TxRx.Ntx,[]);
     symb = mean(:,1) + 1j*mean(:,2);
     
     [~,pos_out1] = min(abs(symb - TxRx.Constellations),[],2);
     pos_out = pos_out1*ones(1,Iter);
     num = sum(pos_out1 ~= pos_in);
end