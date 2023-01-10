
function pos_out = ADAM(TxRx,y,H,N0,HTH,Iter,Constellations,pos_in)
%   AMP detector in Massive MIMO
 sym = TxRx.Constellations(pos_in);
 x = [real(sym.');imag(sym.')];
 m = 2*TxRx.Nrx;
 n = 2*TxRx.Ntx;
%  beta = n/m;
%
 beta1 = 0.9;
 beta2 = 0.99;
 eta = 3e-3;%学习率
    r = y;
    s = zeros(n,1);
    v = zeros(size(HTH,1),size(x,2));
%     tau_old = beta/N0*sigmas2;   %initial estimation variance
    for t=1:Iter
    disp("iter number is "+ Iter);
    g = HTH*x - 2* H.'*y;
    m = beta1*m+(1-beta1)*g;
    v = beta2*v+(1-beta2)*g.^2;
    mhat = m/(1-beta1);
    vhat = v/(1-beta2);
    x = x - eta*vhat/sqrt(mhat);
    beta1 = beta1 *exp(-1);
    beta2 = beta2 *exp(-1);
    end
    mean = reshape(s,TxRx.Ntx,[]);
    symb = mean(:,1) + 1j*mean(:,2);
     
    [~,pos_out1] = min(abs(symb - TxRx.Constellations),[],2);
    pos_out = pos_out1*ones(1,Iter);
    num = sum(pos_out1 ~= pos_in);
end