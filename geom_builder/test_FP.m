
% ts = linspace(-1,1,1000); ts = ts(1:end-1);
% [rsub,dsub,d2sub] = geom_eval(ts,coefs);

ncheb =7;
tsub = linspace(-1,1,1000); tsub = tsub(1:end-1);
% A = [cos(pi*ts.'*(0:ncheb)), sin(pi*ts.'*(1:ncheb))];
% 
% cs = A\d2sub.';
% 
% A1 = [pi*ts.', sin(pi*ts.'*(1:ncheb))./(1:ncheb), -cos(pi*ts.'*(1:ncheb))./(1:ncheb)]/pi;
% A1 = [0*ts.', sin(pi*ts.'*(1:ncheb))./(1:ncheb), -cos(pi*ts.'*(1:ncheb))./(1:ncheb)]/pi;
% 
% cs(1) = 1*pi;
% 
% % cs(1) = 3.2;
% 
% d2cheb = (A*cs).';
% dcheb = (A1*cs).';
% 
% 
% A2 = [pi*ts.', -cos(pi*ts.'*(1:ncheb))./(1:ncheb).^2, -sin(pi*ts.'*(1:ncheb))./(1:ncheb).^2]/pi/pi;
% 
% rcheb = (A2*cs).';
% rcheb = rcheb - mean(rcheb,2) + mean(rsub,2);
% 
% A3 = [0*ts.', -cos(pi*ts.'*(1:ncheb))./(1:ncheb).^2, -sin(pi*ts.'*(1:ncheb))./(1:ncheb).^2]/pi/pi;
% cs2 = A3\(rsub - ts.*[1;0]).';
% cs2(1,1) = pi;
% cs2(1,2) = 0;
% % rcheb = (A2*cs2).';
% % rcheb = rcheb - mean(rcheb,2) + mean(rsub,2);
% 


rs = [(-1:0.3:-0.1).' 0*(-1:0.3:-0.1).'+1; -0.1 0.8; -0.5 0.6; -2/3 -1/2; 0 -1; 2/3 -1/2; 2/3 1/2; 1/2 0.8; (0.5:0.1:1).' (0.5:0.1:1).'*0+1;].';
rs = rs(:,1:end-1);
ts = linspace(-1,1,length(rs)+1); ts = ts(1:end-1);
% ts = [0,cumsum(vecnorm(diff(rs,[],2)))];
% ts = 2*ts/ts(end)-1;
% rs = rs(:,1:end-1);
% ts = ts(1:end-1);

% tsub = ts;
A4 = [1+0*ts.',cos(pi*ts.'*(1:ncheb)), sin(pi*ts.'*(1:ncheb))];

cs4 = A4\(rs - ts.*[1;0]).';

% cs4 = A4\(rsub - ts.*[1;0]).';
rcheb = (A4*cs4).' + ts.*[1;0];
% rcheb = FP2(tsub,cs4);
% rcheb = rcheb - mean(rcheb,2) + mean(rs,2);

% figure(1);clf
% plot(d2sub(1,:), d2sub(2,:),'.')
% hold on
% plot(d2cheb(1,:), d2cheb(2,:),'.')
% hold off
% 
% figure(2);clf
% plot(dsub(1,:), dsub(2,:),'.')
% hold on
% plot(dcheb(1,:), dcheb(2,:),'.')
% hold off
% axis equal

figure(3);clf
plot(rs(1,:), rs(2,:),'x')
hold on
plot(rcheb(1,:), rcheb(2,:),'.')
hold off
axis equal

% tend = [-1,1;];
% Aend = [pi*tend.', -cos(pi*tend.'*(1:ncheb))./(1:ncheb).^2, -sin(pi*tend.'*(1:ncheb))./(1:ncheb).^2]/pi/pi;
% 
% rend = (Aend*cs).';
% d = diff(rend(1,:));


rend = FP2([-1,1],cs4);
d = diff(rend(1,:));


figure(4);clf
plot(rcheb(1,:), rcheb(2,:),'.')
hold on
plot(rcheb(1,:)+d, rcheb(2,:),'.')
plot(rcheb(1,:)+2*d, rcheb(2,:),'.')
hold off


cparams = []; cparams.ta = -1; cparams.tb = 1;
% chnkr = chunkerfunc(@(t) FP(t,cs),cparams);
chnkr = chunkerfunc(@(t) FP2(t,cs4),cparams);


figure(5);clf
plot(chnkr,'.')
hold on
plot(chnkr + [d;0],'.')
hold off

chnkr.npt

function r = FP(ts, cs)
ncheb = (size(cs,1)-1)/2;
A = [pi*ts(:), -cos(pi*ts(:)*(1:ncheb))./(1:ncheb).^2, -sin(pi*ts(:)*(1:ncheb))./(1:ncheb).^2]/pi/pi;

r = (A*cs).';

end

% function r = FP2(ts, cs)
% ncheb = (size(cs,1))/2;
% A = [cos(pi*ts(:)*(1:ncheb)), sin(pi*ts(:)*(1:ncheb))];
% 
% r = (A*cs).' + ts.*[1;0];
% 
% end


function [r,d,d2] = FP2(ts, cs)
ncheb = (size(cs,1))/2;
r  = ([1+0*ts(:),cos(pi*ts(:)*(1:ncheb)), sin(pi*ts(:)*(1:ncheb))]*cs).' + ts(:).'.*[1;0];
d  = ([0*ts(:),-sin(pi*ts(:)*(1:ncheb)), cos(pi*ts(:)*(1:ncheb))]*cs).'+ [1;0];
d2 = ([0*ts(:),-cos(pi*ts(:)*(1:ncheb)), -sin(pi*ts(:)*(1:ncheb))]*cs).';

end

