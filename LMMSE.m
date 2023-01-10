function [pos_out] = LMMSE(TxRx,N0,HTH,HTy,Iter,Es)

  x_es = (HTH + N0/Es*eye(size(HTH)))\HTy;
  
  mean = reshape(x_es,TxRx.Ntx,[]);
  symb = mean(:,1) + 1j*mean(:,2);
  
  [~,pos_out1] = min(abs(symb - TxRx.Constellations),[],2);
  pos_out = pos_out1*ones(1,Iter);
  
return