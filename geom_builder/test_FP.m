ncheb =7;
tsub = linspace(-1,1,1000); tsub = tsub(1:end-1);


rs = [(-1:0.3:-0.1).' 0*(-1:0.3:-0.1).'+1; -0.1 0.8; -0.5 0.6; -2/3 -1/2; 0 -1; 2/3 -1/2; 2/3 1/2; 1/2 0.8; (0.5:0.1:1).' (0.5:0.1:1).'*0+1;].';
rs = rs(:,1:end-1);
ts = linspace(-1,1,length(rs)+1); ts = ts(1:end-1);

% tsub = ts;
A4 = [1+0*ts.',cos(pi*ts.'*(1:ncheb)), sin(pi*ts.'*(1:ncheb))];

cs4 = A4\(rs - ts.*[1;0]).';
rcheb = (A4*cs4).' + ts.*[1;0];


figure(3);clf
plot(rs(1,:), rs(2,:),'x')
hold on
plot(rcheb(1,:), rcheb(2,:),'.')
hold off
axis equal


rend = FP2([-1,1],cs4);
d = diff(rend(1,:));

figure(4);clf
plot(rcheb(1,:), rcheb(2,:),'.')
hold on
plot(rcheb(1,:)+d, rcheb(2,:),'.')
plot(rcheb(1,:)+2*d, rcheb(2,:),'.')
hold off


cparams = []; cparams.ta = -1; cparams.tb = 1;
chnkr = chunkerfunc(@(t) FP2(t,cs4),cparams);


figure(5);clf
plot(chnkr,'.')
hold on
plot(chnkr + [d;0],'.')
hold off

chnkr.npt


function [r,d,d2] = FP2(ts, cs)
ncheb = (size(cs,1))/2;
r  = ([1+0*ts(:),cos(pi*ts(:)*(1:ncheb)), sin(pi*ts(:)*(1:ncheb))]*cs).' + ts(:).'.*[1;0];
d  = ([0*ts(:),-sin(pi*ts(:)*(1:ncheb)), cos(pi*ts(:)*(1:ncheb))]*cs).'+ [1;0];
d2 = ([0*ts(:),-cos(pi*ts(:)*(1:ncheb)), -sin(pi*ts(:)*(1:ncheb))]*cs).';

end

