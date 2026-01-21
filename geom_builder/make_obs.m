
figure(5)
a = drawpolygon;
rs = a.Position.';
if mod(length(rs),2) == 0; rs = rs(:,1:end-1); end

% % rs = rs(:,1:end-1);
% 

% chnkr = sort(chnkr);
% rs = chnkr.r(:,:);
% rs = rs(:,1:end-1);
ts = linspace(-1,1,length(rs)+1); ts = ts(1:end-1);

ncheb = floor((size(rs,2)-1)/2);
% ncheb =7;


% tsub = ts;
A4 = [1+0*ts.',cos(pi*ts.'*(1:ncheb)), sin(pi*ts.'*(1:ncheb))];

cs4 = A4\(rs).';
rcheb = (A4*cs4).';

%%

rcheb = FP2(linspace(-1,1,100), cs4);

figure(3);clf
plot(rs(1,:), rs(2,:),'x')
hold on
plot(rcheb(1,:), rcheb(2,:),'-')
hold off
axis equal

%%
rend = FP2([-1,1],cs4);
% d = diff(rend(1,:));
d = 2;

figure(4);clf
plot(rcheb(1,:), rcheb(2,:),'.')
hold on
plot(rcheb(1,:)+d, rcheb(2,:),'.')
plot(rcheb(1,:)+2*d, rcheb(2,:),'.')
hold off


cparams = []; cparams.ta = -1; cparams.tb = 1; cparams.maxchunklen = 1/(pi/d);
chnkr = chunkerfunc(@(t) FP2(t,cs4),cparams);


figure(5);clf
plot(chnkr,'.')
hold on
for i = 1:3
plot(chnkr + i*[d;0],'.')
end
hold off
axis equal
chnkr.npt

cgrph = tochunkgraph(chnkr);

figure(6);clf
plot(cgrph,'.')
hold on
plot(chnkr,'x')
hold off
axis equal

figure(7);clf
plot_regions(cgrph)

% function [r,d,d2] = FP2(ts, cs)
function [r] = FP2(ts, cs)
ncheb = (size(cs,1))/2;
r  = ([1+0*ts(:),cos(pi*ts(:)*(1:ncheb)), sin(pi*ts(:)*(1:ncheb))]*cs).';
d  = ([0*ts(:),-sin(pi*ts(:)*(1:ncheb)).*(1:ncheb)*pi, cos(pi*ts(:)*(1:ncheb)).*(1:ncheb)*pi]*cs).';
d2 = ([0*ts(:),-cos(pi*ts(:)*(1:ncheb)).*(1:ncheb).^2*pi^2, -sin(pi*ts(:)*(1:ncheb)).*(1:ncheb).^2*pi^2]*cs).';

end

