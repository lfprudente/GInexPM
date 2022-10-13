function [phi] = evalphi(U,V,W,gamma,erropt)

if ( erropt == 1 )
    phi = gamma(1) * norm( V - U,'fro' )^2 + gamma(2) * norm( W - V,'fro' )^2 + gamma(3) * norm( W - U,'fro' )^2;
elseif ( erropt == 2 )
    phi = gamma(1) * norm( V - U,'fro' )^2 + gamma(3) * norm( W - U,'fro' )^2;
elseif ( erropt == 3 )
    phi = gamma(2) * norm( W - V,'fro' )^2 + gamma(3) * norm( W - U,'fro' )^2;
elseif ( erropt == 4 )
    phi = gamma(3) * norm( W - U,'fro' )^2;
elseif ( erropt == 5 )
    phi = gamma(1) * norm( V - U,'fro' )^2 + gamma(2) * norm( W - V,'fro' )^2;
elseif ( erropt == 6 )
    phi = gamma(1) * norm( V - U,'fro' )^2;
elseif ( erropt == 7 )
    phi = gamma(2) * norm( W - V,'fro' )^2;
elseif ( erropt == 8 )
    phi = 0.0;
end