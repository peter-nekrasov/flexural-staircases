func = @(targ) sin(targ(1,:) + 0*targ(2,:));
func = @(targ) (targ(1,:) + 0*targ(2,:));

boxl = 1;
norder = 10;

targ = [0.1;0.0];

[val,grad] = der_by_cheb(func, targ, norder, boxl);

val - func(targ)

grad

grad - cos(targ(1,:) + targ(2,:))

% hess(1) + func(targ)
