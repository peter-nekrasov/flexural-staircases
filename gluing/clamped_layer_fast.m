function u = clamped_layer_fast(chnkr_tr,dens,targ,chnkr,ifree,itrdata,zk,nu,kappa,d,ws,sys,sn,l,rhsmat,rhsmat0,layermat,layermat0)
nkappa = length(kappa);

targmod = [];
targmod.r = real([mod(real(targ.r(1,:)+d/2-mean(chnkr.r(1,:))),d)-d/2+mean(chnkr.r(1,:));targ.r(2,:)]);
% targmod = targ;
nshift = round(real(targ.r(1,:)-targmod.r(1,:))/d);
targmod.r = targ.r(:,:) - nshift.*[d;0];


%%
wts = chnkr_tr.wts(:).'; wts = repmat(wts,4,1);
rhs = -rhsmat*(rhsmat0*(dens .* wts(:)));

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end,:) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end,:))*ws(i);
end

%%

if itrdata == 0
    if ~isempty(layermat)
        u = layermat*(layermat0*sol);
    else
    u = zeros(size(targ.r(:,:),2), 1);
    ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval',kappa,d,sn,[],[],l,0,nu);
    ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval',nu);
    
    wts = repmat(chnkr.wts(:).',2,1);
    
    start1 = tic;
    nbatch = ceil(2e5/chnkr.npt);
    ntout = size(targ.r(:,:),2);
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
        u(iuse,:) = gevalmat*sol;
    end
    t2 = toc(start1);
    fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)
    end
    if ifree
        evalmat = chunkerkernevalmat(chnkr_tr,@(s,t) direct_layer(s,t,zk),targ);
        u = u + evalmat*dens;
    end
else
    if ~isempty(layermat)
        u = layermat*(layermat0*sol);
    else
        ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval_trx',kappa,d,sn,[],[],l,1,nu);
        wts = repmat(chnkr.wts(:).',2,1);
    
        gevalmat = ikern(chnkr,targmod).* wts(:).';
        gevalmat = reshape(gevalmat,nkappa, 4,size(targmod.r,2), []);
        gevalmat = exp(1i*kappa(:).*reshape(nshift,1,1,[])*d) .* gevalmat;
        gevalmat = reshape(gevalmat,nkappa, 4*size(targmod.r,2), []);
        gevalmat = reshape(permute(gevalmat, [2,1,3]), 4*size(targmod.r,2),[]);
        u = gevalmat*sol;
    end
    if ifree
        error('self layer potential not implemented')
    end
end
end