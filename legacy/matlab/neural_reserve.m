function [Pnewseeds, Pnewleaves, Phototropism] = neural_reserve(nb_leaves,vol_relative, Nreserve)
Nhidden	= 3;
Nin 	= 2;
Nout 	= 3;

toto = tan( (Nreserve - 0.5) * pi * 0.99 );
M1 = reshape( toto(1:Nhidden*Nin), Nhidden, Nin);
M2 = reshape( toto(Nhidden*Nin+1:Nhidden*(Nin+Nout)+Nout), Nout, Nhidden+1);
M1(1,1)=0;
M1(Nhidden,Nin)=0;
Xin(1,1) = 0.01 * nb_leaves;
Xin(2,1) = vol_relative;

Zo = M1*Xin;
Zp(1:Nhidden,1) = tanh( 5 * Zo ) / Nhidden;
Zp(Nhidden+1,1) = 1 / Nhidden;
F = M2*Zp;

Pnewleaves   = min( max(0.0,F(1)+2.0), 4.0 ) / 4.0;
Pnewseeds    = min( max(0.0,F(2)+2.0), 4.0 ) / 4.0;
Phototropism = min( max(0.0,F(3)+2.0), 4.0 ) / 4.0;
sum_p = Pnewseeds + Pnewleaves;
if (sum_p > 1)
    Pnewleaves = Pnewleaves / sum_p;
    Pnewseeds  = Pnewseeds  / sum_p;
end
end
