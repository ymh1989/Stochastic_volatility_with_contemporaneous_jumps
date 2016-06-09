clc; clf; clear;
format long
T = 1.0; K = 100; Smax = 300; Vmax = 1;
Nt = 5; dt = T/Nt;
Nx = 30; Nv = 40;
%%% For uniform step
x = linspace(0,Smax,Nx+1); v = linspace(0,Vmax,Nv+1);
% %%%% For random step
% dummy_h_x = rand(1,Nx);
% dummy_h_x(2) = dummy_h_x(1); dummy_h_x(end) = dummy_h_x(end-1);
% h_x = dummy_h_x / sum(dummy_h_x)*Smax;
% dummy_h_v = rand(1,Nv);
% dummy_h_v(2) = dummy_h_v(1); dummy_h_v(end) = dummy_h_v(end-1);
% h_v = dummy_h_v / sum(dummy_h_v)*Vmax;
% x = zeros(1,Nx+1); v = zeros(1,Nv+1);
% for i = 2:Nx+1
%     x(i) = x(i-1) + h_x(i-1); 
% end
% for i = 2:Nv+1
%     v(i) = v(i-1) + h_v(i-1);
% end
%%%%%%%%%%%%%%
% Init
hx = [diff(x) x(end)-x(end-1)]; hv = [diff(v) v(end)-v(end-1)];
main_x(1:Nx-1) = 0; sub_x(1:Nx-1) = 0; sup_x(1:Nx-1) = 0;
main_y(1:Nv-1) = 0; sub_y(1:Nv-1) = 0; sup_y(1:Nv-1) = 0;
fx = zeros(1,Nx-1); fy = zeros(1,Nv-1);
n = 50; zs = linspace(eps, 3, n+1); zv = linspace(0, 3, n+1);

r = 0.03; % riskless rate
q = 0.0; % dividend
sig = 0.3; % vol of vol
rho = -0.5; % correlation between dw1 and dw2
kap = 2; % rate of reversion to the mean level
theta = 0.04; % mean level of var
del = 0.4;
lam = 0.2; % intensity of Poisson process
nu = 0.2; % jump size
gam = -0.5;
rhoj = -0.5;
xi = exp(gam+del^2/2)*(1/(1-nu*rhoj))-1;
p = zeros(n+1); % pdf
for i = 1:n+1
    for j = 1:n+1
        p(i,j) = (1/(sqrt(2*pi)*zs(i)*del*nu)) ...
            * exp(-(zv(j)/(nu)) - ((log(zs(i))-gam-rhoj*zv(j))^2/(2*del^2)));
    end
end

u0 = max(x-K,0)'*ones(1,Nv+1);
[xx0, vv0] = meshgrid(x,v);
%  mesh(xx0,vv0,u0')

x = [x zs(end)*x(end)]; v = [v zv(end)+v(end)];
hx = [hx(1:end-1) x(end)-x(end-1)]; hv = [hv(1:end-1) v(end)-v(end-1)];

% Alternative for interpl2
idx_x(1:Nx+1,1:n+1) = -1;  idxmn_x(1:Nx+1,1:n+1) = -1;
idx_v(1:Nv+1,1:n+1) = -1;  idxmn_v(1:Nv+1,1:n+1) = -1;
for i = 2:Nx+1
    y_x = x(i) * (zs);
    for j = 1:n
        idx_x(i,j) = min(find(x > y_x(j))) - 1; % min(find(x > y(j))) : x_{ij}
        idxmn_x(i,j) = (x(idx_x(i,j)+1) - y_x(j))/hx(idx_x(i,j)); % 1-m_ij
    end
end
for i = 2:Nv+1
    y_v = v(i) + (zv);
    for j = 1:n
        idx_v(i,j) = min(find(v > y_v(j))) - 1; % min(find(x > y(j))) : x_{ij}
        idxmn_v(i,j) = (v(idx_v(i,j)+1) - y_v(j))/hv(idx_v(i,j)); % 1-m_ij
    end
end
% break
% Test for interpl2
for i = 2:Nx+1
    y_x(i-1,:) = x(i) * (zs);
end
for i = 2:Nv+1
    y_v(i-1,:) = v(i) + (zv);
end
% initial condition
u = max(x-K,0)'*ones(1,Nv+2);
% Test
% u = 1/900*(x.^2)'*ones(1,Nv+2);
[xx, vv] = meshgrid(x,v);
% u = interp2(xx0,vv0,u0',xx,vv,'linear');
% mesh(xx,vv,u')

% Test for pdf
% [Xq,Yq] = meshgrid(y_x(Nx-1,:), y_v(Nv-1,:));
% testsimp2D(p',eps, 3, 0, 3, n, n) % integral
% uq = interp2(xx,vv,u',Xq,Yq); % interplation
% mesh(Xq,Yq,uq);
% mesh(Xq,Yq,p');
% mesh(Xq,Yq,uq.*p');
% axis tight
% break;
ou = u;
ui = u;
tic
% break
for l=1:Nt
    ou(1,2:Nv+1) = 0.0;
    ou(end,2:Nv+1) = ou(end-1,2:Nv+1) + hx(end)/hx(end-1)*(ou(end-1,2:Nv+1)-ou(end-2,2:Nv+1));
    ou(:,1) = ou(:,2) + hv(1)/hv(2)*(ou(:,2)-ou(:,3));
    ou(:,end) = ou(:,end-1) + hv(end)/hv(end-1)*(ou(:,end-1)-ou(:,end-2));
    
    for j=2:Nv+1
        for i = 2:Nx+1
            sub_x(i-1) = (-x(i)^2*v(j)...
                +(r-q-lam*xi)*x(i)*hx(i)) / (hx(i-1)*(hx(i-1)+hx(i)));
            main_x(i-1) = 1/dt + r/2 ...
                + (x(i)^2*v(j)-(r-q-lam*xi)*x(i)*(hx(i)-hx(i-1))) / (hx(i-1)*hx(i));
            sup_x(i-1) = (-x(i)^2*v(j)...
                -(r-q-lam*xi)*x(i)*hx(i-1)) / (hx(i)*(hx(i-1)+hx(i)));
            
            % Jump
            [Xq,Yq] = meshgrid(y_x(i-1,:), y_v(j-1,:));
            uq = interp2(xx,vv,ou',Xq,Yq);
            integ = testsimp2D(uq.*p',eps, 3, 0, 3, n, n);
            fx(i-1) = (1/dt-0.5*lam)*ou(i,j)...
                +0.5*rho*sig*x(i)*v(j)...
                *(ou(i-1,j-1)+ou(i+1,j+1)...
                -ou(i+1,j-1)-ou(i-1,j+1))...
                /((hv(j)+hv(j-1))*(hx(i)+hx(i-1))) ...
                + lam*integ;
        end
        %%% Linear boundary condition
        sub_x(end) = sub_x(end) - (hx(end)/hx(end-1))*sup_x(end);
        main_x(end) = main_x(end) + (1+(hx(end)/hx(end-1)))*sup_x(end);
        
        u(2:Nx+1,j)=Thomas(sub_x,main_x,sup_x,fx);
    end
    %%% Linear boundary condition
    u(1,2:Nv+1) = 0.0;
    u(end,2:Nv+1) = u(end-1,2:Nv+1) + hx(end)/hx(end-1)*(u(end-1,2:Nv+1)-u(end-2,2:Nv+1));
    u(:,1) = u(:,2) + hv(1)/hv(2)*(u(:,2)-u(:,3));
    u(:,end) = u(:,end-1) + hv(end)/hv(end-1)*(u(:,end-1)-u(:,end-2));
    for i=2:Nx+1
        for j = 2:Nv+1
            sub_y(j-1) = (-sig^2*v(j)...
                +kap*(theta-v(j))*hv(j))/(hv(j-1)*(hv(j-1)+hv(j)));
            main_y(j-1) = 1/dt + r/2 + ...
                (sig^2*v(j)-kap*(theta-v(j))*(hv(j)-hv(j-1)))/(hv(j-1)*hv(j));
            sup_y(j-1) = (-sig^2*v(j)...
                -kap*(theta-v(j))*hv(j-1))/(hv(j)*(hv(j-1)+hv(j)));
            
            % Jump
            [Xq,Yq] = meshgrid(y_x(i-1,:), y_v(j-1,:));
            uq = interp2(xx,vv,ou',Xq,Yq);
            integ = testsimp2D(uq.*p',eps, 3, 0, 3, n, n);
            
            fy(j-1) = (1/dt-0.5*lam)*u(i,j)...
                +0.5*rho*sig*x(i)*v(j)...
                *(u(i-1,j-1)+u(i+1,j+1)...
                -u(i+1,j-1)-u(i-1,j+1))...
                /((hv(j)+hv(j-1))*(hx(i)+hx(i-1))) ...
                + lam*integ;
        end
        % Linear boundary condition
        main_y(1) = main_y(1) + (1+(hv(1)/hv(2)))*sub_y(1);
        sup_y(1) = sup_y(1)-(hv(1)/hv(2))*sub_y(1);
        sub_y(end) = sub_y(end) - (hv(end)/hv(end-1))*sup_y(end);
        main_y(end) = main_y(end) + (1+(hv(end)/hv(end-1)))*sup_y(end);
        
        ou(i,2:Nv+1) = Thomas(sub_y,main_y,sup_y,fy);
    end    
end
%%% Linear boundary condition
ou(1,2:Nv+1) = 0.0;
ou(end,2:Nv+1) = ou(end-1,2:Nv+1) + hx(end)/hx(end-1)*(ou(end-1,2:Nv+1)-ou(end-2,2:Nv+1));
ou(:,1) = ou(:,2) + hv(1)/hv(2)*(ou(:,2)-ou(:,3));
ou(:,end) = ou(:,end-1) + hv(end)/hv(end-1)*(ou(:,end-1)-ou(:,end-2));
toc
mesh(xx(1:end-1,1:end-1),vv(1:end-1,1:end-1),ou(1:end-1,1:end-1)')
% mesh(xx,vv,ou');
axis tight
view(-45, 30)
title('Numerical solution of SVCJ')
xlabel('S');
ylabel('vol');
zlabel('call price')