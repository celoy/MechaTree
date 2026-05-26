function [Safety] = neural_branch(nb_leaves, max_stress, Nbranch)
Nhidden	= 3;
Nin 	= 2;
Nout 	= 1;

toto = tan( (Nbranch - 0.5) * pi * 0.99 );
M1 = reshape( toto(1:Nhidden*Nin), Nhidden, Nin);
M2 = reshape( toto(Nhidden*Nin+1:Nhidden*(Nin+Nout)+Nout), Nout, Nhidden+1);
M1(1,1)=0;
M1(Nhidden,Nin)=0;
Xin(1,1) = 0.01 * nb_leaves;
Xin(2,1) = max_stress;

Zo = M1*Xin;
Zp(1:Nhidden,1) = tanh( 5 * Zo ) / Nhidden;
Zp(Nhidden+1,1) = 1 / Nhidden;
F = M2*Zp;
Safety   = max(0.0,F(1)+1.0);
end
