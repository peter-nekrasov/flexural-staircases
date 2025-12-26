function [r,d,d2] = geom_eval(tsub,coefs)

as = coefs{1};
bs = coefs{2};
cs = coefs{3}; 
ds = coefs{4};

ts = linspace(-1,1,length(as)+1);

inds = sum(tsub(:) >= ts,2);
inds(inds == length(as)+1) = length(as);

dx = (tsub(:).' - ts(inds)).' ;

r = as(inds,:) + bs(inds,:).*dx+ cs(inds,:).*dx.^2 + ds(inds,:).*dx.^3;
d = bs(inds,:) + 2*cs(inds,:).*dx + 3*ds(inds,:).*dx.^2;
d2 = 2*cs(inds,:) + 6*ds(inds,:).*dx;

% npol = length(cs);
% 
% tsub = linspace(-1,1,length(tsub));
% 
% % A_r = lege.pols(ts,npol-1);
% % A_r = squeeze(A_r).';
% 
% A_r = cos(acos(tsub(:)).*(0:npol-1));

% A_r = [ts(:).^(0:1) cos(d*(1:npol).*ts(:))];
% A_d = [(0:npol-1).*ts(:).^(-1:npol-2) -2*d*(1:npol).*sin(2*d*(1:npol).*ts(:))];
% A_d2 = [(-1:npol-2).*(0:npol-1).*ts(:).^(-2:npol-3) -(2*d*(1:npol)).^2.*cos(2*d*(1:npol).*ts(:))];

% r = A_r*cs;
% d = A_d*cs;
% d2 = A_d2*cs;

r = r.';
d = d.';
d2 = d2.';

end
