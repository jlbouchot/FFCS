function [prox_TV, prox_TVstar] = proximal_TV(z,sigma)

prox_TV = sign2(z).*(abs2(z)-sigma).*(abs2(z)>=sigma) ;
prox_TVstar = z-sigma*prox_F(1/sigma, z/sigma);

return