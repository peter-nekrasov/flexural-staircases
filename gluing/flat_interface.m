function [f,fd,fdd] = flat_interface(t, a, b, t0, t1)
    
    phi   = @(t,u,v,z) u*(t-z).*erfc(u*(t-z))*v - exp(-u^2*(t-z).^2)/sqrt(pi)*v;
    phid  = @(t,u,v,z) u*erfc(u*(t-z))*v;
    phidd = @(t,u,v,z) -u*u*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    f = zeros([2,size(t)]);
    fd = zeros([2,size(t)]);
    fdd = zeros([2,size(t)]);
    
    f(1,:) = t + 1i*(phi(t,a,b,t0) - phi(t,-a,b,t1)); 
    fd(1,:)= 1 + 1i*(phid(t,a,b,t0) - phid(t,-a,b,t1));
    fdd(1,:) = 1i*(phidd(t,a,b,t0) - phidd(t,-a,b,t1));
    
    f(2,:) = 0;
    fd(2,:) = 0;
    fdd(2,:) = 0;
        
end

