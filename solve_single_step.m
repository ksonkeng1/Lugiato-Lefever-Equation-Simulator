function uf = solve_single_step(uf,ELOf,Linv_Ff,NL,h)

ut = ifft(ELOf.*(uf + Linv_Ff) - Linv_Ff);
uf2 = fft(ut.*NL(ut,h));
uf = ELOf.*(uf2 + Linv_Ff)-Linv_Ff;

end
