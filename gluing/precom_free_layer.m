function [rhsmat,layermat] = precom_free_layer(chnkr_tr,targ,chnkr,itrdata,zk,nu,kappa,d,sn,l,s0_l,sn_l)

ising = 1;
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs_trx',kappa,d,sn,s0_l,sn_l,l,ising,nu);
rhsmat = bskern(chnkr_tr,chnkr);

nkappa = length(kappa);

targmod = [];
targmod.r = real([mod(real(targ.r(1,:)+d/2),d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round(real(targ.r(1,:)-targmod.r(1,:))/d);
targmod.r = targ.r(:,:) - nshift.*[d;0];


if ~itrdata
    ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',kappa,d,sn,s0_l,sn_l,l,0,nu);
    ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval',nu);
    
    wts = repmat(chnkr.wts(:).',3,1);
    
    ntout = size(targ.r(:,:),2);
    gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targmod);
    layermat = ikern(chnkr,targmod).* wts(:).';
    
    layermat = reshape(layermat,nkappa, ntout, []);
    layermat = layermat + reshape(gevalmat_0,1,ntout, []);
    layermat = exp(1i*kappa(:).*nshift*d) .* layermat;
    layermat = reshape(permute(layermat, [2,1,3]), ntout,[]);
else
    ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval_trx',kappa,d,sn,s0_l,sn_l,l,1,nu);
    wts = repmat(chnkr.wts(:).',3,1);

    layermat = ikern(chnkr,targmod).* wts(:).';
    layermat = reshape(layermat,nkappa, 4,size(targmod.r,2), []);
    layermat = exp(1i*kappa(:).*reshape(nshift,1,1,[])*d) .* layermat;
    layermat = reshape(layermat,nkappa, 4*size(targmod.r,2), []);
    layermat = reshape(permute(layermat, [2,1,3]), 4*size(targmod.r,2),[]);
end
end