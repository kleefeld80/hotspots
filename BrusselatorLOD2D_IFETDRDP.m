function [runtime,w_old] = BrusselatorLOD2D_IFETDRDP(dt,steps,ende,fignr,filename,no)

% dt: time step
% steps: number of spatial points in each coordinate direction

%% Model Paramters and initial conditions
A = 1;

% diffusion coefficient
epsln1 = 0.1;

% create nodes
x = linspace(0,1,steps); 
h = abs(x(1)-x(2));
fprintf('dt=%f h=%f\n',dt,h)
y = x;
nnodes = steps^2;
nodes = zeros(nnodes,2);
j = 1;
for k = 1 : steps
        for i = 1:steps
               nodes(j,:) = [x(i) y(k)];
            j = j+1;
        end
end
nb = nnodes; % becuase we are solving a syste of 2 RDE

% discretize time interval
t = 0:dt:ende; tlen = length(t)

% initial condition for u
randn('seed',10);
u_old=randn(length(nodes),1);

% Stacking nodes for evolution
w_old = zeros(nb,1);
w_old = u_old;

%% Block matrix Assembly
C = zeros(1);
C(1,1) = (epsln1*dt)/h^2;
Q = blktridiag(2,-1,-1,steps);
Q(1,2) = -2; Q(steps,steps-1) = -2; I = speye(steps);
A1 = kron(kron(I,Q),C);A2 = kron(kron(Q,I),C); 
I = speye(nb); 
M1 = sparse(I+A1); M2 = sparse(I+(1/3)*A1); M3 = sparse(I+(1/4)*A1);
M11 = sparse(I+A2); M22 = sparse(I+(1/3)*A2); M33 = sparse(I+(1/4)*A2);

%% Time Evolution 
[L1,U1] =lu(M1);
[L2,U2] =lu(M2);
[L3,U3] =lu(M3);
 [L11,U11] =lu(M11); 
 [L22,U22] =lu(M22);
 [L33,U33] =lu(M33);

tic
disp('Halll')
for i = 2:tlen
     F_old = F(w_old);
     w_star = U11\(L11\(w_old + dt*F_old));
     w_star = U1\(L1\w_star);
     F_star = F(w_star);
     a_1 = U2\(L2\w_old);
     b_1 = U3\(L3\w_old);
     c_1 = 9*a_1 - 8*b_1;
     a_2 = U2\(L2\F_old);
     b_2 = U3\(L3\F_old);
     c_2 = 9*a_2-8*b_2;
     d_1 = U22\(L22\(9*c_1+2*dt*c_2 + dt*F_star));
     d_2 = U33\(L33\(-8*c_1-(3/2)*dt*c_2-(dt/2)*F_star));
     w_old = d_1+d_2;              
end
disp('Wer')
%w_old
 u_soln = w_old; 
  U = reshape(u_soln,steps,steps);
  
  terrain = [0 166 0
        22 174 0
        45 182 0
        71 190 0
        99 198 0
        128 206 0
        160 214 0
        194 222 0
        230 230 0
        231 208 26
        232 192 52
        234 182 78
        235 177 105
        237 178 131
        238 185 159
        239 198 186
        241 217 214
        242 242 242];
    terrain = terrain/255;
    
  figure(fignr)
       hold on
       ticks = 12; % Schriftgroesse der Achsenunterteilungen
       set(gca,'fontsize',ticks,'fontweight','bold','LineWidth',2);
        %Z = Z/Z(min(find(abs(Z(:)) == max(abs(Z(:))))));
          contourf(x,x,U,no)
          cb=colorbar;
        cb.Ticks = linspace(0,1,5);
        colormap(terrain)
        colorbar('off')
        %h.LineStyle ='none';
        xticks([0,0.2,0.4,0.6,0.8,1])
        yticks([0,0.2,0.4,0.6,0.8,1])
        size = 12; % Schriftgroesse der Achsenbeschriftung
        xlabel('x','FontSize',size,'fontweight','bold');
        ylabel('y','FontSize',size,'fontweight','bold');
        axis equal
        box on
        hold off
        print(filename,'-depsc2')
  
runtime = toc;

function Fr = F(U)
 Fr = zeros(nb,1);
 u = U(1:2:nb-1); v = U(2:2:nb);
 f1 = 0;
 f2 = 0;
 %Fr(1:2:nb-1) = f1; Fr(2:2:nb) = f2;
end

end
