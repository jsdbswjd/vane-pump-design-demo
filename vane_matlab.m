% Patent-informed vane pump C1 reconstruction, MATLAB R2018b+. No toolbox.
% Assumptions: OD53.5, wall2, rotor radius19 mm, 12 zero-thickness vanes.
clear;clc;K=2;alpha=10;width=14;rpm=3000;target=13;
z=calc(K,alpha,width,rpm,1440);
fprintf('Displacement %.6f cc/rev, ideal flow %.6f L/min\n',z.disp,z.flow);
fprintf('Rmin %.6f mm, max pressure angle %.3f deg, feasible %d\n',z.rmin,z.pa,z.valid);
figure('Color','w');
subplot(1,3,1);plot(z.r.*cos(z.phi),z.r.*sin(z.phi),'LineWidth',2);hold on;
plot(19*cos(z.phi),19*sin(z.phi));plot(26.75*cos(z.phi),26.75*sin(z.phi));axis equal;grid on;
xlabel('x [mm]');ylabel('y [mm]');title('C1 reconstructed cam');
subplot(1,3,2);plot(rad2deg(z.phi(1:end-1)),z.v);grid on;
xlabel('Rotor angle [deg]');ylabel('One chamber volume [cc]');
ks=linspace(1.5,3,61);ds=nan(size(ks));
for i=1:numel(ks),s=calc(ks(i),alpha,width,rpm,720);if s.valid,ds(i)=s.disp;end,end
subplot(1,3,3);plot(ks,ds);hold on;yline(target,'--');grid on;
xlabel('Profile K');ylabel('Displacement [cc/rev]');title('Feasible K sweep');
% Grid minimum width for requested displacement (not global optimum).
best=[NaN NaN Inf];
for k=1.5:0.05:3
 for a=4:0.5:15
  s=calc(k,a,1,rpm,720);b=max(8,target/s.disp);
  if s.valid && b<=18 && b<best(3),best=[k a b];end
 end
end
if isfinite(best(3)),fprintf('Grid best K %.2f, angle %.1f deg, width %.4f mm\n',best);else,fprintf('No feasible grid design\n');end
function z=calc(k,alpha,b,rpm,N)
 a=deg2rad(alpha);R=24.75/(k+2);t=2*a;
 T=[R*(t-sin(t))-pi*R,R*(1-cos(t)+k)];direction=[sin(a),cos(a)];
 F=T-dot(T,direction)*direction;rmin=norm(F);beta=atan2(F(2),-F(1));
 u=linspace(0,beta,100)';arc=[-rmin*cos(u),rmin*sin(u)];
 u=linspace(0,1,100)';line=F+u.*(T-F);
 t=linspace(2*a,pi,400)';cycl=[R*(t-sin(t))-pi*R,R*(1-cos(t)+k)];
 pts=[arc;line(2:end,:);cycl(2:end,:)];aa=atan2(pts(:,2),-pts(:,1));rr=sqrt(sum(pts.^2,2));
 phi=linspace(0,2*pi,N+1);be=atan2(abs(sin(phi)),abs(cos(phi)));r=interp1(aa,rr,be,'linear');
 area=.5*(r.^2-19^2);seg=(area(1:end-1)+area(2:end))/2*(2*pi/N)*b;
 cs=[0,cumsum([seg seg])];idx=1:N;v=(cs(idx+N/12)-cs(idx))/1000;
 displacement=24*(max(v)-min(v));pa=max(rad2deg(atan(abs(gradient(r,2*pi/N)./r))));
 z=struct('phi',phi,'r',r,'v',v,'disp',displacement,'flow',displacement*rpm/1000,'rmin',rmin,'pa',pa,'valid',rmin-19>=.5 && pa<=30);
end
