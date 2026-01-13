function u = free_layer(chnkr_tr,dens,targ,chnkr,ifree,itrdata,zk,nu,kappa,d,ws,sys,sn,l,H,s0_l,sn_l)
nkappa = length(kappa);

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);


%%
ising = 1;
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs_trx',kappa,d,sn,s0_l,sn_l,l,ising,nu);

wts = chnkr_tr.wts(:).';
wts = repmat(wts,4,1);
dens = dens .* wts(:);

rhs = -bskern(chnkr_tr,chnkr)*dens;

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end,:) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end,:))*ws(i);
end
sol = reshape(sol,nkappa,[],size(sol,2));

dens_comb = zeros(nkappa,3*chnkr.npt,size(rhs,2));
dens_comb(:,1:3:end,:) = sol(:,1:2:end,:);
tmpsol = permute(sol,[2,3,1]);
dh = pagemtimes(H,tmpsol(1:2:end,:,:));
dh = permute(dh,[3,1,2]);
dens_comb(:,2:3:end,:) = dh;
dens_comb(:,3:3:end,:) = sol(:,2:2:end,:);
dens_comb = reshape(dens_comb,[],size(rhs,2));

%%
if itrdata == 0
    u = zeros(size(targ.r,2), size(src.r,2));
    
    ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',kappa,d,sn,s0_l,sn_l,l,0,nu);
    ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval',nu);
    
    wts = repmat(chnkr.wts(:).',3,1);
    
    start1 = tic;
    nbatch = ceil(2e5/chnkr.npt);
    ntout = size(targ.r,2);
    for i = 1:ceil(ntout/nbatch)
        iuse = ((i-1)*nbatch+1):min(ntout,i*nbatch);
        targi = []; targi.r = targmod.r(:,iuse);
        nti = length(iuse);
        
        gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targi);
        gevalmat = ikern(chnkr,targi).* wts(:).';
        
        gevalmat = reshape(gevalmat,nkappa, nti, []);
        gevalmat = gevalmat + reshape(gevalmat_0,1,nti, []);
        gevalmat = exp(1i*kappa(:).*nshift(iuse)*d) .* gevalmat;
        % gevalmat = exp(1i*kappa(:).*nshift(iout)*d) .* gevalmat;
        
        gevalmat = reshape(permute(gevalmat, [2,1,3]), nti,[]);
        u(iuse,:) = gevalmat*dens_comb;
    end
    t2 = toc(start1);
    fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

    if ifree
        u = u + ...
            chunkerkerneval(chnkr_tr,@(s,t) direct_layer(s,t,zk),dens,targ);
    end
else
    ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval_trx',kappa,d,sn,s0_l,sn_l,l,1,nu);
    wts = repmat(chnkr.wts(:).',3,1);

    gevalmat = ikern(chnkr,targmod).* wts(:).';
    gevalmat = reshape(gevalmat,nkappa, 4,size(targmod.r,2), []);
    gevalmat = exp(1i*kappa(:).*reshape(nshift,1,1,[])*d) .* gevalmat;
    gevalmat = reshape(gevalmat,nkappa, 4*size(targmod.r,2), []);
    gevalmat = reshape(permute(gevalmat, [2,1,3]), 4*size(targmod.r,2),[]);
    u = gevalmat*dens_comb;

    if ifree
        error('self layer potential not implemented')
    end
end
end