function [val,grad,hess,third,fourth] = green(src,targ,zk,kappa,d,sn,l,ising)
%CHNK.FLEX2DQUAS.GREEN evaluate the quasiperiodic flexural Green's function
% for the given sources and targets
%
% Input:
%   zk - wavenumber
%   kappa - quasiperiodic parameters
%   d - period
%   sn - precomputed lattice sum integrals 
%       (see chnk.helm2dquas.latticecoefs)
%   l - number of periodic copies computed explicitly
%   ising - if set to 0, only include the periodic copies. If set to 1,
%       include the free-space part
%
% see also CHNK.FLEX2DQUAS.KERN
[~,nsrc] = size(src);
[~,ntarg] = size(targ);

xs = repmat(src(1,:),ntarg,1);
ys = repmat(src(2,:),ntarg,1);

xt = repmat(targ(1,:).',1,nsrc);
yt = repmat(targ(2,:).',1,nsrc);

rx = xt-xs;
ry = yt-ys;


rx = rx(:);
ry = ry(:);


nx = fix(rx/d);
rx = rx - nx*d;


rx2 = rx.*rx;
ry2 = ry.*ry;

r2 = rx2+ry2;

r = sqrt(r2);

npt = size(r,1);

ythresh = 2*d/2;
iclose = abs(ry) < ythresh;
ifar = ~iclose;

rxfar = rx(ifar);
ryfar = ry(ifar);

rxclose = rx(iclose);
ryclose = ry(iclose);
rclose = r(iclose);

nxclose = nx(iclose);
nptclose = size(rxclose, 1);

nkappa = length(kappa);

val = zeros(nkappa,npt,1);
if nargout > 1
grad = zeros(nkappa,npt,2);
end
if nargout > 2
hess = zeros(nkappa,npt,3);
end
if nargout > 3
third = zeros(nkappa,npt,4);
end
if nargout > 4
fourth = zeros(nkappa,npt,5);
end

tol = 1e-10;
Lbd = sqrt((log(tol))^2/real(ythresh)^2 + real(zk)^2);

rxfar = rxfar.';
ryfar = ryfar.';
if ~isempty(ryfar)
M = ceil(Lbd*d/(2*pi));
ms = reshape((-M:M),1,1,[]);
xi_m = kappa(:) + 2*pi/d*ms;

% beta = sqrt((xi_m.^2-zk^2));
betah = sqrt(1i*(xi_m-zk)).*sqrt(-1i*(xi_m+zk));
betak = sqrt(1i*(xi_m-1i*zk)).*sqrt(-1i*(xi_m+1i*zk));

fhath = exp(-betah.*sqrt(ryfar.^2) + 1i*xi_m.*rxfar)./(2*betah);
fhatk = exp(-betak.*sqrt(ryfar.^2) + 1i*xi_m.*rxfar)./(2*betak);

xifac = 1i*xi_m;
bfach = -betah.*(sqrt(ryfar.^2)./ryfar);
bfack = -betak.*(sqrt(ryfar.^2)./ryfar);

val(:,ifar,:) = sum(fhath - fhatk,3)/(d);
if nargout > 1
grad(:,ifar,1) = sum(xifac.*(fhath-fhatk),3)/d;
grad(:,ifar,2) = sum(bfach.*fhath-bfack.*fhatk,3)/d;
end

if nargout >2
hess(:,ifar,1) = sum(xifac.^2.*(fhath-fhatk),3)/d;
hess(:,ifar,2) = sum(xifac.*(bfach.*fhath-bfack.*fhatk),3)/d;
hess(:,ifar,3) = sum(bfach.^2.*fhath-bfack.^2.*fhatk,3)/d;
end

if nargout > 3
third(:,ifar,1) = sum(xifac.^3.*(fhath-fhatk),3)/d;
third(:,ifar,2) = sum(xifac.^2.*(bfach.*fhath-bfack.*fhatk),3)/d;
third(:,ifar,3) = sum(xifac.*(bfach.^2.*fhath-bfack.^2.*fhatk),3)/d;
third(:,ifar,4) = sum(bfach.^3.*fhath-bfack.^3.*fhatk,3)/d;
end

if nargout > 4
fourth(:,ifar,1) = sum(xifac.^4.*(fhath-fhatk),3)/d;
fourth(:,ifar,2) = sum(xifac.^3.*(bfach.*fhath-bfack.*fhatk),3)/d;
fourth(:,ifar,3) = sum(xifac.^2.*(bfach.^2.*fhath-bfack.^2.*fhatk),3)/d;
fourth(:,ifar,4) = sum(xifac.*(bfach.^3.*fhath-bfack.^3.*fhatk),3)/d;
fourth(:,ifar,5) = sum(bfach.^4.*fhath-bfack.^4.*fhatk,3)/d;
end

end

alpha = (exp(1i*kappa(:)*d));

val_near= zeros(nkappa,nptclose);
grad_near = zeros(nkappa,nptclose,2);
hess_near = zeros(nkappa,nptclose,3);
third_near = zeros(nkappa,nptclose,4);
fourth_near = zeros(nkappa,nptclose,5);
ls = -l:l;
if ~isempty(rxclose)
    for i = ls
        if ising == 1
            iuse = true(nptclose,1);
        else
            iuse = nxclose ~= -i;
        end

        rxi = rxclose - i*d;
        if nargout>4
        [vali,gradi,hessi,thirdi,fourthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;   
        third_near(:,iuse,:) = third_near(:,iuse,:) + thirdi(:,iuse,:).*alpha.^i;   
        fourth_near(:,iuse,:) = fourth_near(:,iuse,:) + fourthi(:,iuse,:).*alpha.^i;   
        elseif nargout>3
        [vali,gradi,hessi,thirdi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;   
        third_near(:,iuse,:) = third_near(:,iuse,:) + thirdi(:,iuse,:).*alpha.^i;        
        elseif nargout>2
        [vali,gradi,hessi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;
        elseif nargout > 1
        [vali,gradi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        else
        vali = chnk.flex2d.hkdiffgreen(zk,[0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        end
    end

    N = size(sn,2)-1;
    ns = (0:N);
    ns_use = (0:N+4);
    Js = zeros(length(rclose),N+5);
    Is = zeros(length(rclose),N+5);

    if length(rclose) < N+5
        for i = 1:length(rclose)
        Js(i,:) = besselj(ns_use,zk*rclose(i));
        Is(i,:) = besselj(ns_use,1i*zk*rclose(i));
        end
    else
        for i = 1:length(ns_use)
        Js(:,i) = besselj(ns_use(i),zk*rclose);
        Is(:,i) = besselj(ns_use(i),1i*zk*rclose);
        end
    end
    % t1 = tic;
    eip = (rxclose+1i*ryclose)./rclose;
    eipn = reshape(eip.^ns,1,[], N+1);
    cs = (eipn+1./eipn)/2;
    
    Js = reshape(Js,1,[],N+5);
    Is = reshape(Is,1,[],N+5);
    snj = reshape(sn(:,:,1),nkappa, 1, N+1);
    sni = reshape(sn(:,:,2),nkappa, 1, N+1);
    
    tmpj = reshape(Js(:,:,2:end-4).*cs(:,:,2:end),[],N);
    tmpi = reshape(-Is(:,:,2:end-4).*cs(:,:,2:end),[],N);    
    val_far = 0.25*1i*(Js(:,:,1).*snj(:,:,1)-Is(:,:,1).*sni(:,:,1)) + 0.5*1i*snj(:,2:end)*tmpj.' + 0.5*1i*sni(:,2:end)*tmpi.';
    val(:,iclose) = val_near+val_far;
    
    if nargout >1
        DJs = cat(3,-Js(:,:,2),.5*(Js(:,:,1:end-5)-Js(:,:,3:end-3)))*zk;
        DIs = cat(3,-Is(:,:,2),.5*(Is(:,:,1:end-5)-Is(:,:,3:end-3)))*1i*zk;
        ss = (eipn-1./eipn)/2i;
        
        tmpj = reshape(DJs(:,:,2:end).*cs(:,:,2:end),[],N);
        tmpi = reshape(-DIs(:,:,2:end).*cs(:,:,2:end),[],N);

        grad_far_p = 0.25*1i*(DJs(:,:,1).*snj(:,:,1)-DIs(:,:,1).*sni(:,:,1)) + 0.5*1i*snj(:,2:end)*tmpj.' + 0.5*1i*sni(:,2:end)*tmpi.';
        
        tmpj = reshape((Js(:,:,2:end-4)).*ss(:,:,2:end),[],N)./rclose;
        tmpi = reshape((-Is(:,:,2:end-4)).*ss(:,:,2:end),[],N)./rclose;

        grad_far_t = (0.5*1i*((-reshape((1:N),1,[]).*sni(:,2:end))*tmpi.'))+(0.5*1i*((-reshape((1:N),1,[]).*snj(:,2:end))*tmpj.'));
        
        grad_far = cat(3,cs(:,:,2).*grad_far_p - ss(:,:,2).*grad_far_t, ss(:,:,2).*grad_far_p + cs(:,:,2).*grad_far_t);
        grad(:,iclose,:) = grad_near + grad_far; 
    end
    if nargout > 2
        DDJs = cat(3,.5*(Js(:,:,3)-Js(:,:,1)),.25*(Js(:,:,4)-3*Js(:,:,2)),.25*(Js(:,:,1:end-6)-2*Js(:,:,3:end-4)+Js(:,:,5:end-2)))*zk^2;
        DDIs = cat(3,.5*(Is(:,:,3)-Is(:,:,1)),.25*(Is(:,:,4)-3*Is(:,:,2)),.25*(Is(:,:,1:end-6)-2*Is(:,:,3:end-4)+Is(:,:,5:end-2)))*(1i*zk)^2;
        rclose = rclose.';
        rxclose = rxclose.';
        ryclose = ryclose.';
        ns = reshape(ns,1,1,[]);

        tmp_nj = rclose.^(-4).*(-ns.*ryclose.*(Js(:,:,1:end-4)).*(ns.*ryclose.*cs+2*rxclose.*ss)+ ...
            rclose.*ryclose.*(ryclose.*cs + 2*ns.*rxclose.*ss).*(DJs)+ ...
            rclose.^2.*rxclose.^2.*cs.*(DDJs));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-4).*(-ns.*ryclose.*(-Is(:,:,1:end-4)).*(ns.*ryclose.*cs+2*rxclose.*ss)+ ...
            rclose.*ryclose.*(ryclose.*cs + 2*ns.*rxclose.*ss).*(-DIs)+ ...
            rclose.^2.*rxclose.^2.*cs.*(-DDIs));
        tmp_ni = reshape(tmp_ni,[],N+1);
        hess_far_xx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = ns.*(Js(:,:,1:end-4)).*(ns.*rxclose.*ryclose.*cs+(rxclose.^2-ryclose.^2).*ss).*rclose.^(-4) ...
            +rclose.^(-3).*(-(rxclose.*ryclose.*cs+ns.*(rxclose.^2-ryclose.^2).*ss).*(DJs)+...
            rclose.*rxclose.*ryclose.*cs.*(DDJs));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = ns.*(-Is(:,:,1:end-4)).*(ns.*rxclose.*ryclose.*cs+(rxclose.^2-ryclose.^2).*ss).*rclose.^(-4) ...
            +rclose.^(-3).*(-(rxclose.*ryclose.*cs+ns.*(rxclose.^2-ryclose.^2).*ss).*(-DIs)+...
            rclose.*rxclose.*ryclose.*cs.*(-DDIs));
        tmp_ni = reshape(tmp_ni,[],N+1);
        hess_far_xy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-4).*(-ns.*rxclose.*(Js(:,:,1:end-4)).*(ns.*rxclose.*cs-2*ryclose.*ss)+ ...
            rclose.*rxclose.*(rxclose.*cs - 2*ns.*ryclose.*ss).*(DJs)+ ...
            rclose.^2.*ryclose.^2.*cs.*(DDJs));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-4).*(-ns.*rxclose.*(-Is(:,:,1:end-4)).*(ns.*rxclose.*cs-2*ryclose.*ss)+ ...
            rclose.*rxclose.*(rxclose.*cs - 2*ns.*ryclose.*ss).*(-DIs)+ ...
            rclose.^2.*ryclose.^2.*cs.*(-DDIs));
        tmp_ni = reshape(tmp_ni,[],N+1);
        hess_far_yy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        hess_far = cat(3,hess_far_xx, hess_far_xy, hess_far_yy);

        hess(:,iclose,:) = hess_near + hess_far;
    end
    if nargout > 3
        D3Js = cat(3,.75*Js(:,:,2) - .25*Js(:,:,4), ...
            .25*(2*Js(:,:,3)-.5*Js(:,:,5)-1.5*Js(:,:,1)), ...
            .25*(-2*Js(:,:,2)+1.5*Js(:,:,4)-.5*Js(:,:,6)), ...
            .25*(.5*Js(:,:,1:end-7)-1.5*Js(:,:,3:end-5)+1.5*Js(:,:,5:end-3)-.5*Js(:,:,7:end-1))) *zk^3;
        D3Is = cat(3,.75*Is(:,:,2) - .25*Is(:,:,4), ...
            .25*(2*Is(:,:,3)-.5*Is(:,:,5)-1.5*Is(:,:,1)), ...
            .25*(-2*Is(:,:,2)+1.5*Is(:,:,4)-.5*Is(:,:,6)), ...
            .25*(.5*Is(:,:,1:end-7)-1.5*Is(:,:,3:end-5)+1.5*Is(:,:,5:end-3)-.5*Is(:,:,7:end-1))) *(1i*zk)^3;

        tmp_nj = rclose.^(-6).*(ns.*ryclose.*(Js(:,:,1:end-4)).*(6*ns.*rxclose.*ryclose.*cs+(6*rxclose.^2 - (2+ns.^2).*ryclose.^2).*ss)+ ...
            rclose.*(-3.*ryclose.*((1+ns.^2).*rxclose.*ryclose.*cs + ns.*(2*rxclose.^2-ryclose.^2).*ss)).*(DJs)+ ...
            rclose.^2.*3.*rxclose.*ryclose.*(ryclose.*cs+ns.*rxclose.*ss).*(DDJs) + ...
            rclose.^3.*rxclose.^3.*cs.*(D3Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*ryclose.*(-Is(:,:,1:end-4)).*(6*ns.*rxclose.*ryclose.*cs+(6*rxclose.^2 - (2+ns.^2).*ryclose.^2).*ss)+ ...
            rclose.*(-3.*ryclose.*((1+ns.^2).*rxclose.*ryclose.*cs + ns.*(2*rxclose.^2-ryclose.^2).*ss)).*(-DIs)+ ...
            rclose.^2.*3.*rxclose.*ryclose.*(ryclose.*cs+ns.*rxclose.*ss).*(-DDIs) + ...
            rclose.^3.*rxclose.^3.*cs.*(-D3Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_xxx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';


        tmp_nj = rclose.^(-6).*(ns.*(Js(:,:,1:end-4)).*(2*ns.*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + rxclose.*(-2*rxclose.^2+(6+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + ns.*rxclose.*(2*rxclose.^2-7*ryclose.^2).*ss).*(DJs) + ...
            -rclose.^2.*(-ryclose.*(-2*rxclose.^2 + ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 2*ryclose.^2).*ss).*(DDJs) + ...
            rclose.^3.*(rxclose.^2.*ryclose).*cs.*(D3Js)); 
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*(-Is(:,:,1:end-4)).*(2*ns.*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + rxclose.*(-2*rxclose.^2+(6+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*ryclose.*(-2*rxclose.^2+ryclose.^2).*cs + ns.*rxclose.*(2*rxclose.^2-7*ryclose.^2).*ss).*(-DIs) + ...
            -rclose.^2.*(-ryclose.*(-2*rxclose.^2 + ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 2*ryclose.^2).*ss).*(-DDIs) + ...
            rclose.^3.*(rxclose.^2.*ryclose).*cs.*(-D3Is)); 
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_xxy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-6).*(ns.*(Js(:,:,1:end-4)).*(2*ns.*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs + ryclose.*(2*ryclose.^2-(6+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs - ns.*ryclose.*(2*ryclose.^2-7*rxclose.^2).*ss).*(DJs) + ...
            rclose.^2.*(rxclose.*(-2*ryclose.^2 + rxclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 2*rxclose.^2).*ss).*(DDJs) + ...
            rclose.^3.*(ryclose.^2.*rxclose).*cs.*(D3Js)); 
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*(-Is(:,:,1:end-4)).*(2*ns.*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs + ryclose.*(2*ryclose.^2-(6+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-(1+ns.^2).*rxclose.*(-2*ryclose.^2+rxclose.^2).*cs - ns.*ryclose.*(2*ryclose.^2-7*rxclose.^2).*ss).*(-DIs) + ...
            rclose.^2.*(rxclose.*(-2*ryclose.^2 + rxclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 2*rxclose.^2).*ss).*(-DDIs) + ...
            rclose.^3.*(ryclose.^2.*rxclose).*cs.*(-D3Is)); 
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_xyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-6).*(ns.*rxclose.*(Js(:,:,1:end-4)).*(6*ns.*ryclose.*rxclose.*cs+(-6*ryclose.^2 + (2+ns.^2).*rxclose.^2).*ss)+ ...
            rclose.*(-3.*rxclose.*((1+ns.^2).*ryclose.*rxclose.*cs + ns.*(-2*ryclose.^2+rxclose.^2).*ss)).*(DJs)+ ...
            rclose.^2.*3.*ryclose.*rxclose.*(rxclose.*cs-ns.*ryclose.*ss).*(DDJs) + ...
            rclose.^3.*ryclose.^3.*cs.*(D3Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-6).*(ns.*rxclose.*(-Is(:,:,1:end-4)).*(6*ns.*ryclose.*rxclose.*cs+(-6*ryclose.^2 + (2+ns.^2).*rxclose.^2).*ss)+ ...
            rclose.*(-3.*rxclose.*((1+ns.^2).*ryclose.*rxclose.*cs + ns.*(-2*ryclose.^2+rxclose.^2).*ss)).*(-DIs)+ ...
            rclose.^2.*3.*ryclose.*rxclose.*(rxclose.*cs-ns.*ryclose.*ss).*(-DDIs) + ...
            rclose.^3.*ryclose.^3.*cs.*(-D3Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        third_far_yyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        third_far = cat(3,third_far_xxx,third_far_xxy,third_far_xyy,third_far_yyy);

        third(:,iclose,:) = third_near + third_far;
    end
    if nargout > 4
        D4Js = cat(3,.75/2*(Js(:,:,1)-Js(:,:,3)) - .25/2*(Js(:,:,3)-Js(:,:,5)), ...
            .25*((Js(:,:,2)-Js(:,:,4))-.5/2*(Js(:,:,4)-Js(:,:,6))+1.5*Js(:,:,2)), ...
            .25*(-(Js(:,:,1)-Js(:,:,3))+1.5/2*(Js(:,:,3)-Js(:,:,5))-.5/2*(Js(:,:,5)-Js(:,:,7))), ...
            .25*(-.5*Js(:,:,2)-1.5/2*(Js(:,:,2)-Js(:,:,4))+1.5/2*(Js(:,:,4)-Js(:,:,6))-.5/2*(Js(:,:,6)-Js(:,:,8))), ...
            .25/2*(.5*(Js(:,:,1:end-8)-Js(:,:,3:end-6))-1.5*(Js(:,:,3:end-6)-Js(:,:,5:end-4))+1.5*(Js(:,:,5:end-4)-Js(:,:,7:end-2))-.5*(Js(:,:,7:end-2)-Js(:,:,9:end)))) *zk^4;
        D4Is = cat(3,.75/2*(Is(:,:,1)-Is(:,:,3)) - .25/2*(Is(:,:,3)-Is(:,:,5)), ...
            .25*((Is(:,:,2)-Is(:,:,4))-.5/2*(Is(:,:,4)-Is(:,:,6))+1.5*Is(:,:,2)), ...
            .25*(-(Is(:,:,1)-Is(:,:,3))+1.5/2*(Is(:,:,3)-Is(:,:,5))-.5/2*(Is(:,:,5)-Is(:,:,7))), ...
            .25*(-.5*Is(:,:,2)-1.5/2*(Is(:,:,2)-Is(:,:,4))+1.5/2*(Is(:,:,4)-Is(:,:,6))-.5/2*(Is(:,:,6)-Is(:,:,8))), ...
            .25/2*(.5*(Is(:,:,1:end-8)-Is(:,:,3:end-6))-1.5*(Is(:,:,3:end-6)-Is(:,:,5:end-4))+1.5*(Is(:,:,5:end-4)-Is(:,:,7:end-2))-.5*(Is(:,:,7:end-2)-Is(:,:,9:end))))  *(1i*zk)^4;

        tmp_nj = rclose.^(-8).*(ns.*ryclose.*(Js(:,:,1:end-4)).*(ns.*ryclose.*(-36*rxclose.^2+(8+ns.^2).*ryclose.^2).*cs + 12*rxclose.*(-2*rxclose.^2+(2+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(ryclose.*(-3*(1+2*ns.^2).*ryclose.*(-4*rxclose.^2+ryclose.^2).*cs - 4*ns.*rxclose.*(-6*rxclose.^2 + (8+ns.^2).*ryclose.^2).*ss).*(DJs)) + ...
            rclose.^2.*(3*ryclose.*((-2*(2+ns.^2).*rxclose.^2.*ryclose + ryclose.^3).*cs + 4*ns.*rxclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(DDJs) + ...
            rclose.^3.*rxclose.^2.*(2*ryclose.*(3*ryclose.*cs + 2*ns.*rxclose.*ss)).*(D3Js) + ...
            rclose.^4.*rxclose.^4.*cs.*(D4Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*ryclose.*(-Is(:,:,1:end-4)).*(ns.*ryclose.*(-36*rxclose.^2+(8+ns.^2).*ryclose.^2).*cs + 12*rxclose.*(-2*rxclose.^2+(2+ns.^2).*ryclose.^2).*ss) + ...
            rclose.*(ryclose.*(-3*(1+2*ns.^2).*ryclose.*(-4*rxclose.^2+ryclose.^2).*cs - 4*ns.*rxclose.*(-6*rxclose.^2 + (8+ns.^2).*ryclose.^2).*ss).*(-DIs)) + ...
            rclose.^2.*(3*ryclose.*((-2*(2+ns.^2).*rxclose.^2.*ryclose + ryclose.^3).*cs + 4*ns.*rxclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(-DDIs) + ...
            rclose.^3.*rxclose.^2.*(2*ryclose.*(3*ryclose.*cs + 2*ns.*rxclose.*ss)).*(-D3Is) + ...
            rclose.^4.*rxclose.^4.*cs.*(-D4Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xxxx = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(ns.*(Js(:,:,1:end-4)).*(-ns.*rxclose.*ryclose.*(-18*rxclose.^2 + (26+ns.^2).*ryclose.^2).*cs + 3*(2*rxclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*ryclose.^4).*ss) ...
            - rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(2*rxclose.^2-3*ryclose.^2).*cs + ns.*(6*rxclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*ryclose.^4).*ss).*(DJs) ...
            + 3*rclose.^2.*(rxclose.*ryclose.*((2+ns.^2).*rxclose.^2 - (3+ns.^2).*ryclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(DDJs) - ...
            rxclose.*rclose.^3.*((3*ryclose.*(rxclose.^2-ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 3*ryclose.^2).*ss).*(D3Js) - rxclose.^2.*ryclose.*rclose.*cs.*(D4Js)));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*(-Is(:,:,1:end-4)).*(-ns.*rxclose.*ryclose.*(-18*rxclose.^2 + (26+ns.^2).*ryclose.^2).*cs + 3*(2*rxclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*ryclose.^4).*ss) ...
            - rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(2*rxclose.^2-3*ryclose.^2).*cs + ns.*(6*rxclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*ryclose.^4).*ss).*(-DIs) ...
            + 3*rclose.^2.*(rxclose.*ryclose.*((2+ns.^2).*rxclose.^2 - (3+ns.^2).*ryclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(-DDIs) - ...
            rxclose.*rclose.^3.*((3*ryclose.*(rxclose.^2-ryclose.^2).*cs + ns.*rxclose.*(rxclose.^2 - 3*ryclose.^2).*ss).*(-D3Is) - rxclose.^2.*ryclose.*rclose.*cs.*(-D4Is)));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xxxy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(ns.*(Js(:,:,1:end-4)).*(ns.*(-6*rxclose.^4 + (32+ns.^2).*rxclose.^2.*ryclose.^2 - 6*ryclose.^4).*cs + 6*(4+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss) + ...
            rclose.*((1+2*ns.^2).*(2*rxclose.^4 - 11*rxclose.^2.*ryclose.^2 + 2*ryclose.^4).*cs - 2*ns.*(14+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss).*(DJs) - ...
            rclose.^2.*(((2+ns.^2).*rxclose.^4 - (11+4*ns.^2).*rxclose.^2.*ryclose.^2 + (2+ns.^2).*ryclose.^4).*cs + 12*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*(DDJs) + ...
            rclose.^3.*((rxclose.^4 - 4*rxclose.^2.*ryclose.^2 + ryclose.^4).*cs + 2*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*(D3Js) + ...
            rclose.^4.*(rxclose.^2.*ryclose.^2.*cs).*(D4Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*(-Is(:,:,1:end-4)).*(ns.*(-6*rxclose.^4 + (32+ns.^2).*rxclose.^2.*ryclose.^2 - 6*ryclose.^4).*cs + 6*(4+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss) + ...
            rclose.*((1+2*ns.^2).*(2*rxclose.^4 - 11*rxclose.^2.*ryclose.^2 + 2*ryclose.^4).*cs - 2*ns.*(14+ns.^2).*rxclose.*ryclose.*(rxclose.^2 - ryclose.^2).*ss).*(-DIs) - ...
            rclose.^2.*(((2+ns.^2).*rxclose.^4 - (11+4*ns.^2).*rxclose.^2.*ryclose.^2 + (2+ns.^2).*ryclose.^4).*cs + 12*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*(-DDIs) + ...
            rclose.^3.*((rxclose.^4 - 4*rxclose.^2.*ryclose.^2 + ryclose.^4).*cs + 2*ns.*rxclose.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss).*( - D3Is) + ...
            rclose.^4.*(rxclose.^2.*ryclose.^2.*cs).*(-D4Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xxyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(-ns.*(Js(:,:,1:end-4)).*(ns.*rxclose.*ryclose.*(-18*ryclose.^2 + (26+ns.^2).*rxclose.^2).*cs + 3*(2*ryclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*rxclose.^4).*ss) ...
            + rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(-2*ryclose.^2+3*rxclose.^2).*cs + ns.*(6*ryclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*rxclose.^4).*ss).*(DJs) ...
            - 3*rclose.^2.*(rxclose.*ryclose.*(-(2+ns.^2).*ryclose.^2 + (3+ns.^2).*rxclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(DDJs) + ...
            ryclose.*rclose.^3.*((3*rxclose.*(rxclose.^2-ryclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 3*rxclose.^2).*ss).*(D3Js) + ryclose.^2.*rxclose.*rclose.*cs.*(D4Js)));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(-ns.*(-Is(:,:,1:end-4)).*(ns.*rxclose.*ryclose.*(-18*ryclose.^2 + (26+ns.^2).*rxclose.^2).*cs + 3*(2*ryclose.^4 - 3*(4+ns.^2).*rxclose.^2.*ryclose.^2+(2+ns.^2).*rxclose.^4).*ss) ...
            + rclose.*(3*(1+2*ns.^2).*rxclose.*ryclose.*(-2*ryclose.^2+3*rxclose.^2).*cs + ns.*(6*ryclose.^4-3*(14+ns.^2).*rxclose.^2.*ryclose.^2 + (8+ns.^2).*rxclose.^4).*ss).*(-DIs) ...
            - 3*rclose.^2.*(rxclose.*ryclose.*(-(2+ns.^2).*ryclose.^2 + (3+ns.^2).*rxclose.^2).*cs + ns.*(rxclose.^4 - 6*rxclose.^2.*ryclose.^2 + ryclose.^4).*ss).*(-DDIs) + ...
            ryclose.*rclose.^3.*((3*rxclose.*(rxclose.^2-ryclose.^2).*cs + ns.*ryclose.*(ryclose.^2 - 3*rxclose.^2).*ss).*(-D3Is) + ryclose.^2.*rxclose.*rclose.*cs.*(-D4Is)));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_xyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';

        tmp_nj = rclose.^(-8).*(ns.*rxclose.*(Js(:,:,1:end-4)).*(ns.*rxclose.*(-36*ryclose.^2+(8+ns.^2).*rxclose.^2).*cs - 12*ryclose.*(-2*ryclose.^2+(2+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-rxclose.*(3*(1+2*ns.^2).*rxclose.*(-4*ryclose.^2+rxclose.^2).*cs - 4*ns.*ryclose.*(-6*ryclose.^2 + (8+ns.^2).*rxclose.^2).*ss).*(DJs)) + ...
            rclose.^2.*(3*rxclose.*((-2*(2+ns.^2).*ryclose.^2.*rxclose + rxclose.^3).*cs + 4*ns.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(DDJs) + ...
            rclose.^3.*ryclose.^2.*(2*rxclose.*(3*rxclose.*cs - 2*ns.*ryclose.*ss)).*(D3Js) + ...
            rclose.^4.*ryclose.^4.*cs.*(D4Js));
        tmp_nj = reshape(tmp_nj,[],N+1);
        tmp_ni = rclose.^(-8).*(ns.*rxclose.*(-Is(:,:,1:end-4)).*(ns.*rxclose.*(-36*ryclose.^2+(8+ns.^2).*rxclose.^2).*cs - 12*ryclose.*(-2*ryclose.^2+(2+ns.^2).*rxclose.^2).*ss) + ...
            rclose.*(-rxclose.*(3*(1+2*ns.^2).*rxclose.*(-4*ryclose.^2+rxclose.^2).*cs - 4*ns.*ryclose.*(-6*ryclose.^2 + (8+ns.^2).*rxclose.^2).*ss).*(-DIs)) + ...
            rclose.^2.*(3*rxclose.*((-2*(2+ns.^2).*ryclose.^2.*rxclose + rxclose.^3).*cs + 4*ns.*ryclose.*(-rxclose.^2 + ryclose.^2).*ss)).*(-DDIs) + ...
            rclose.^3.*ryclose.^2.*(2*rxclose.*(3*rxclose.*cs - 2*ns.*ryclose.*ss)).*(-D3Is) + ...
            rclose.^4.*ryclose.^4.*cs.*(-D4Is));
        tmp_ni = reshape(tmp_ni,[],N+1);
        fourth_far_yyyy = 0.25*1i*tmp_nj(:,1).'.*snj(:,1)+.5*1i*snj(:,2:end)*tmp_nj(:,2:end).' + ...
            0.25*1i*tmp_ni(:,1).'.*sni(:,1)+.5*1i*sni(:,2:end)*tmp_ni(:,2:end).';        

        fourth_far = cat(3,fourth_far_xxxx,fourth_far_xxxy,fourth_far_xxyy,fourth_far_xyyy,fourth_far_yyyy);

        fourth(:,iclose,:) = fourth_near + fourth_far;
    end
end

quasi_phase = exp(1i*kappa(:)*nx(:).'*d);


if nargout == 1
    val = quasi_phase.*val;

    if ising == 0
        isub = (abs(nx(:)) > max(ls)) | ifar;

        if any(isub)
        vali = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).'+ nx(isub).'*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        val(:,isub,:) = val(:,isub,:) - vali;
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
elseif nargout == 2
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    
    if ising == 0
        isub = (abs(nx(:)) > max(ls)) | ifar;
    
        if any(isub)
        [vali, gradi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        val(:,isub,:) = val(:,isub,:) - vali;
        grad(:,isub,:) = grad(:,isub,:) - gradi;
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
elseif nargout == 3
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    
    if ising == 0
        isub = (abs(nx(:)) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);

        val(:,isub,:) = val(:,isub,:) - vali;
        grad(:,isub,:) = grad(:,isub,:) - gradi;
        hess(:,isub,:) = hess(:,isub,:) - hessi;
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
elseif nargout == 4
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    third = quasi_phase.*third;
    
    if ising == 0
        isub = (abs(nx(:)) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi,thirdi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);

        val(:,isub,:) = val(:,isub,:) - vali;
        grad(:,isub,:) = grad(:,isub,:) - gradi;
        hess(:,isub,:) = hess(:,isub,:) - hessi;
        third(:,isub,:) = third(:,isub,:) - thirdi;
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
    third = reshape(third,nkappa*ntarg,nsrc,4);
elseif nargout == 5
    val = quasi_phase.*val;
    grad = quasi_phase.*grad;
    hess = quasi_phase.*hess;
    third = quasi_phase.*third;
    fourth = quasi_phase.*fourth;
    
    if ising == 0
        isub = (abs(nx(:)) > max(ls)) | ifar;

        if any(isub)
        [vali, gradi, hessi,thirdi,fourthi] = chnk.flex2d.hkdiffgreen(zk,[0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        fourthi = reshape(fourthi,1,[],5);

        val(:,isub,:) = val(:,isub,:) - vali;
        grad(:,isub,:) = grad(:,isub,:) - gradi;
        hess(:,isub,:) = hess(:,isub,:) - hessi;
        third(:,isub,:) = third(:,isub,:) - thirdi;
        fourth(:,isub,:) = fourth(:,isub,:) - fourthi;
        end
    end

    val = reshape(val,nkappa*ntarg,nsrc);
    grad = reshape(grad,nkappa*ntarg,nsrc,2);
    hess = reshape(hess,nkappa*ntarg,nsrc,3);
    third = reshape(third,nkappa*ntarg,nsrc,4);
    fourth = reshape(fourth,nkappa*ntarg,nsrc,5);
end
end