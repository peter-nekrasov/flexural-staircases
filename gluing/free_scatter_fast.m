function u = free_scatter_fast(src,targ,chnkr,ifree,itrdata,zk,nu,kappa,d,ws,sys,sn,l,H,s0_l,sn_l,layermat,layermat0)
% itrdata = 0, get field, itrdata = 1, get field and first 3 x_1 ders
nkappa = length(kappa);


targmod = [];
targmod.r = real([mod(real(targ.r(1,:)+d/2-mean(chnkr.r(1,:))),d)-d/2+mean(chnkr.r(1,:));targ.r(2,:)]);
% targmod = targ;
nshift = round(real(targ.r(1,:)-targmod.r(1,:))/d);
targmod.r = targ.r(:,:) - nshift.*[d;0];


%%
ising = 1;
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa,d,sn,s0_l,sn_l,l,ising,nu);

rhs = -bskern(src,chnkr);

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

u = layermat*(layermat0*dens_comb);

if itrdata == 0
    if ifree
        u = u + chnk.flex2d.kern(zk,src,targ, 's');
    end
else
    if ifree
        u = u + direct_rhs(src,targ,zk);
    end
end
end