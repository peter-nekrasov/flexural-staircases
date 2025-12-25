function coefs = get_splines(rs)

nspline = length(rs)-1;
ts = linspace(-1,1,nspline+1);
h = ts(2)-ts(1);

% ts = cos(pi*((0:(length(rs)-1))+1/2)/length(rs));
% A = cos(acos(ts(:)).*(0:npol-1));
% A = lege.pols(ts,npol-1);
% A = squeeze(A).';
% A = [A; (1 - (-1).^(1:npol)).*(0:npol-1).*(1:npol)/2];
% cs = A \ rs.'; %; [0 0]];

% vecnorm(rs.' - (A*cs))

c = zeros(3*nspline,1);
r = zeros(1,3*nspline);
r(1) = h; c(1) = h;
r(nspline+1) = h^2;
r(2*nspline+1) = h^3;

A1 = toeplitz(c,r);
A1 = A1(1:nspline,:);

c = zeros(3*nspline,1);
r = zeros(1,3*nspline);
r(1) = 1; c(1) = 1;
r(2) = -1;
r(nspline+1) = 2*h;
r(2*nspline+1) = 3*h^2;

A2 = toeplitz(c,r);
A2 = A2(1:nspline-1,:);

% c = zeros(3*nspline,1);
r = zeros(1,3*nspline);
r(nspline+1) = 2;
r(nspline+2) = -2;
r(2*nspline+1) = 6*h;

A3 = toeplitz(r);
A3 = A3(1:nspline-1,:);

A4 = zeros(2,3*nspline);
A4(1,1) = 1;
A4(2,end-2*nspline) = 1;
A4(2,end-nspline) = 2*h;
A4(2,end) = 3*h^2;

A = [A1; A2; A3; A4];

rhs = zeros(3*nspline,2);
rhs(1:nspline,:) = (rs(:,2:end) - rs(:,1:end-1)).';
rhs(end-1:end,1) = 1;

cfs = A \ rhs;

as = rs(:,1:end-1).';
bs = cfs(1:nspline,:);
cs = cfs(nspline+1:2*nspline,:);
ds = cfs(2*nspline+1:end,:);

coefs = {as,bs,cs,ds};

end