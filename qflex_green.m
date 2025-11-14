function [val, grad, hess, third, fourth] = qflex_green(src,targ,zk,kappa,d,pxys,cs,l,ising)

pgreen = @(s,t) chnk.flex2d.hkdiffgreen(zk,s,t);


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


if ~isempty(ryfar)

    if nargout > 4
    [val(:,ifar,:), grad(:,ifar,:), hess(:,ifar,:), third(:,ifar,:), fourth(:,ifar,:)] =  quasi_flex_dual_sum(rxfar,ryfar,zk,kappa,d);
    elseif nargout > 3
    [val(:,ifar,:), grad(:,ifar,:), hess(:,ifar,:), third(:,ifar,:)] =  quasi_flex_dual_sum(rxfar,ryfar,zk,kappa,d);
    elseif nargout > 2
    [val(:,ifar,:), grad(:,ifar,:), hess(:,ifar,:)] =  quasi_flex_dual_sum(rxfar,ryfar,zk,kappa,d);
    elseif nargout > 1
    [val(:,ifar,:), grad(:,ifar,:)] =  quasi_flex_dual_sum(rxfar,ryfar,zk,kappa,d);
    else
    val(:,ifar,:) =  quasi_flex_dual_sum(rxfar,ryfar,zk,kappa,d);
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
        [vali,gradi,hessi,thirdi,fourthi] = pgreen([0;0],[rxi.';ryclose.']);
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
        [vali,gradi,hessi,thirdi] = pgreen([0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        thirdi = reshape(thirdi,1,[],4);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;
        third_near(:,iuse,:) = third_near(:,iuse,:) + thirdi(:,iuse,:).*alpha.^i;
        elseif nargout>2
        [vali,gradi,hessi] = pgreen([0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        hessi = reshape(hessi,1,[],3);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        hess_near(:,iuse,:) = hess_near(:,iuse,:) + hessi(:,iuse,:).*alpha.^i;
        elseif nargout > 1
        [vali,gradi] = pgreen([0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        gradi = reshape(gradi,1,[],2);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        grad_near(:,iuse,:) = grad_near(:,iuse,:) + gradi(:,iuse,:).*alpha.^i;
        else
        vali = pgreen([0;0],[rxi.';ryclose.']);
        vali = reshape(vali,1,[],1);
        val_near(:,iuse) = val_near(:,iuse) + vali(:,iuse).*alpha.^i;
        end
    end

    if nargout > 4
        [val_far, grad_far, hess_far, third_far, fourth_far] = pgreen(pxys,[rxclose.';ryclose.']);
        val_far = reshape(val_far * cs, nptclose, nkappa, 1);
        grad_far = reshape(pagemtimes(grad_far, cs), nptclose, nkappa, 2);
        hess_far = reshape(pagemtimes(hess_far, cs), nptclose, nkappa, 3);
        third_far = reshape(pagemtimes(third_far, cs), nptclose, nkappa, 4);
        fourth_far = reshape(pagemtimes(fourth_far, cs), nptclose, nkappa, 5);

        val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
        grad(:,iclose,:) = grad_near + permute(grad_far, [2,1,3]); 
        hess(:,iclose,:) = hess_near + permute(hess_far, [2,1,3]);
        third(:,iclose,:) = third_near + permute(third_far, [2,1,3]);
        fourth(:,iclose,:) = fourth_near + permute(fourth_far, [2,1,3]);
    elseif nargout > 3
        [val_far, grad_far, hess_far, third_far] = pgreen(pxys,[rxclose.';ryclose.']);
        val_far = reshape(val_far * cs, nptclose, nkappa, 1);
        grad_far = reshape(pagemtimes(grad_far, cs), nptclose, nkappa, 2);
        hess_far = reshape(pagemtimes(hess_far, cs), nptclose, nkappa, 3);
        third_far = reshape(pagemtimes(third_far, cs), nptclose, nkappa, 4);

        val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
        grad(:,iclose,:) = grad_near + permute(grad_far, [2,1,3]); 
        hess(:,iclose,:) = hess_near + permute(hess_far, [2,1,3]);
        third(:,iclose,:) = third_near + permute(third_far, [2,1,3]);
    elseif nargout > 2
        [val_far, grad_far, hess_far] = pgreen(pxys,[rxclose.';ryclose.']);
        val_far = reshape(val_far * cs, nptclose, nkappa, 1);
        grad_far = reshape(pagemtimes(grad_far, cs), nptclose, nkappa, 2);
        hess_far = reshape(pagemtimes(hess_far, cs), nptclose, nkappa, 3);

        val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
        grad(:,iclose,:) = grad_near + permute(grad_far, [2,1,3]); 
        hess(:,iclose,:) = hess_near + permute(hess_far, [2,1,3]);
    elseif nargout >1
        [val_far, grad_far] = pgreen(pxys,[rxclose.';ryclose.']);
        val_far = reshape(val_far * cs, nptclose, nkappa, 1);
        grad_far = reshape(pagemtimes(grad_far, cs), nptclose, nkappa, 2);
        
        val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
        grad(:,iclose,:) = grad_near + permute(grad_far, [2,1,3]);
    else
        val_far = pgreen(pxys,[rxclose.';ryclose.']);
        val_far = reshape(val_far * cs, nptclose, nkappa, 1);
        
        val(:,iclose,:) = val_near + permute(val_far, [2,1,3]); 
    end
end

quasi_phase = exp(1i*kappa(:)*nx(:).'*d);


if nargout == 1
    val = quasi_phase.*val;

    if ising == 0
        isub = (abs(nx(:)) > max(ls)) | ifar;

        if any(isub)
        vali = pgreen([0;0],[rx(isub).'+ nx(isub).'*d;ry(isub).']);
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
        [vali, gradi] = pgreen([0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
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
        [vali, gradi, hessi] = pgreen([0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
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
        [vali, gradi, hessi,thirdi] = pgreen([0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
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
        [vali, gradi, hessi,thirdi,fourthi] = pgreen([0;0],[rx(isub).' + nx(isub).'*d;ry(isub).']);
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
