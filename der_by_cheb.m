function [val,grad,hess,third,fourth,fifth] = der_by_cheb(func, targ, norder, boxl)

ss = cos((2*(1:norder)-1)/2/norder*pi);

[nordersx,nordersy] = meshgrid(0:norder-1);
norders = [nordersx(:).';nordersy(:).'];

Tcheb = cos(nordersx(:).' .*acos(ss(1+nordersx(:)).')).*cos(nordersy(:).' .*acos(ss(1+nordersy(:))).');

xs = boxl*ss;
[X,Y] = meshgrid(xs);
targs = [X(:).'; Y(:).'];

val0 = func(targs);
coefs = Tcheb\val0.';

targ0 = targ/boxl;
Teval = cos(nordersx(:).' .*acos(targ0(1,:))).*cos(nordersy(:).' .*acos(targ0(2,:)));
val = Teval*coefs;

% csx = cos(nordersx(:).' .*acos(targ0(1,:)));
% csy = cos(nordersy(:).' .*acos(targ0(2,:)));
% ssx = sin(nordersx(:).' .*acos(targ0(1,:)));
% ssy = sin(nordersy(:).' .*acos(targ0(2,:)));

dTx = (nordersx(:).'>0).*2.*(nordersx(:).'+1).*Teval;
dTy = (nordersy(:).'>0).*2.*(nordersy(:).'+1).*Teval;


% THAT ISN'T THE DERIVATIVE
if nargout >1
grad = zeros(2,length(val));
grad(1,:) = dTx*coefs;
grad(2,:) = dTy*coefs;
grad = grad/boxl;
end

if nargout >2
hess = zeros(3,length(val));
% hess(1,:) = (-nordersx(:).'.^2.*csx.*csy)*coefs;
% hess(2,:) = (nordersx(:).'.*nordersy(:).'.*ssx.*ssy)*coefs;
% hess(3,:) = (-nordersy(:).'.^2.*csx.*csy)*coefs;
hess = hess/boxl;
end

if nargout >3
third = zeros(3,length(val));
% third(1,:) = (nordersx(:).'.^3.*ssx.*csy)*coefs;
% third(2,:) = (nordersx(:).'.^2.*nordersy(:).'.*csx.*ssy)*coefs;
% third(3,:) = (nordersx(:).'.*nordersy(:).'.^2.*ssx.*csy)*coefs;
% third(4,:) = (-nordersy(:).'.^3.*csx.*ssy)*coefs;
third = third/boxl;
end

if nargout >4
fourth = zeros(3,length(val));
end

if nargout >5
fifth = zeros(3,length(val));
end


end