function u = clamped_scatter(src,targ,chnkr,ifree,itrdata,zk,nu,kappa,d,ws,sys,sn,l)
% itrdata = 0, get field, itrdata = 1, get field and first 3 x1 ders
nkappa = length(kappa);


targmod = [];
targmod.r = real([mod(real(targ.r(1,:)+d/2),d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round(real(targ.r(1,:)-targmod.r(1,:))/d);

targmod.r = targ.r(:,:) - nshift.*[d;0];


%%
ising = 1;
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_bcs',kappa,d,sn,[],[],l,ising,nu);

rhs = -bskern(src,chnkr);

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end,:) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end,:))*ws(i);
end
% sol = reshape(sol,nkappa,[],size(sol,2));

%%

if itrdata == 0
    u = zeros(size(targ.r(:,:),2), size(src.r(:,:),2));
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
    
    if ifree
        u = u + chnk.flex2d.kern(zk,src,targ, 's');
    end
else
    ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval_trx',kappa,d,sn,[],[],l,1,nu);
    wts = repmat(chnkr.wts(:).',2,1);

    gevalmat = ikern(chnkr,targmod).* wts(:).';
    gevalmat = reshape(gevalmat,nkappa, 4,size(targmod.r,2), []);
    gevalmat = exp(1i*kappa(:).*reshape(nshift,1,1,[])*d) .* gevalmat;
    gevalmat = reshape(gevalmat,nkappa, 4*size(targmod.r,2), []);
    gevalmat = reshape(permute(gevalmat, [2,1,3]), 4*size(targmod.r,2),[]);
    u = gevalmat*sol;

    if ifree
        [val,grad,hess,third] = chnk.flex2d.hkdiffgreen(zk,src.r,targ.r);  
        submat = zeros(4*size(targ.r(:,:),2),size(src.r(:,:),2));
        submat(1:4:end,:) = 1/(2*zk^2).*val;
        submat(2:4:end,:) = 1/(2*zk^2).*grad(:,:,1);
        submat(3:4:end,:) = 1/(2*zk^2).*hess(:,:,1);
        submat(4:4:end,:) = 1/(2*zk^2).*third(:,:,1);
        u = u + submat;
    end
end
end