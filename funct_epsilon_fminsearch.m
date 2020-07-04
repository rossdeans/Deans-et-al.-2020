%function to solve for TLP
%now done so that it can be solved using fminsearch not fsolve
function f=funct_epsilon_fminsearch(Ca,g_star,dw,mu_K_mean,mu_chi_mean,mu_piv_mean,psi_s,epsilon)
    f=@F;
    
    function y=F(x)
        y=((mu_piv_mean/mu_chi_mean)*(1+x/epsilon)*(x+psi_s)^2-(1+(mu_K_mean*dw*(1+x/epsilon)/mu_chi_mean)^0.5-psi_s/epsilon)*(((1.6*(1+x/epsilon)*(Ca-g_star)*(x+psi_s))/mu_chi_mean)^0.5-1.6*(1+(mu_K_mean*dw*(1+x/epsilon)/mu_chi_mean)^0.5+x/epsilon))).^2;
    end
end