function HelmDL2alphaBump1Trick(ns)
    global alpha
  global l
  global m
  global N
  global tol_rank
  global mp
  global R
  global faces
  global n
  global vi
  global vhi
  global wavenumber
  
  n=ns;
    alpha=(1-sqrt(3/5))/2;
    faces=n/2;
    fprintf('Number of faces: %d\n',faces);
    fprintf('Number of collocation nodes: %d\n',3*faces);
    
    % parameters
  N = 24; 
  l=20;
  tol_rank = 10^(-4);
  R = 0.5;
  mp = 0.8;
  wavenumber=1;
    
    vi=geometrie3(40);
    
    % scale the points
    vi=vi*0.25;
    n=size(vi,2)/2-1
    faces=n/2;
    
    vhi=geometrie3(20);
    nnn=length(vhi);
    % scale the points
    vhi=vhi*0.25;
    outer=vhi(:,1:nnn/2);
    inner=vhi(:,nnn/2+1:end);
    %reso=100;
        
    %x=linspace(-1.6,1.6,reso);
    %y=linspace(-1.6,1.6,reso);
    %[xx,yy]=meshgrid(x,y);
    %X=load('complexBump1')
    %plotthething2('complexBump1zoom.eps',xx(50:100,40:60),yy(50:100,40:60),real(X.z(50:100,40:60)),outer,inner,12,[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2],[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
    %pause
    %plotthething5('geocomplex.eps',outer,inner,10,[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2],[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
    %pause   
    
    indi1=[1:2:n,n+2:2:2*n];
    indi2=[2:2:n,n+3:2:2*n+1];
    indi3=[3:2:n+1,n+4:2:2*(n+1)];
    
    vhi=vi;
    k=1;
    L=1;
    for i=1:2*faces
      x1=vi(1,indi1(i));
      y1=vi(2,indi1(i));
      x2=vi(1,indi2(i));
      y2=vi(2,indi2(i));
      x3=vi(1,indi3(i));
      y3=vi(2,indi3(i));
      s=alpha;
      vhi(1,L)=l1(s)*x1+l2(s)*x2+l3(s)*x3;
      vhi(2,L)=l1(s)*y1+l2(s)*y2+l3(s)*y3; 
      s=0.5;
      vhi(1,L+1)=l1(s)*x1+l2(s)*x2+l3(s)*x3;
      vhi(2,L+1)=l1(s)*y1+l2(s)*y2+l3(s)*y3; 
      s=1-alpha;
      vhi(1,L+2)=l1(s)*x1+l2(s)*x2+l3(s)*x3;
      vhi(2,L+2)=l1(s)*y1+l2(s)*y2+l3(s)*y3; 
      k=k+2;
      L=L+3;
    end
    
    figure(4)
    hold on
    plot(vi(1,:),vi(2,:));
    plot(vhi(1,:),vhi(2,:),'x');
    hold off
    %pause
    trala=1;
    if trala
    m=2*3*faces; 
    [eigenwerte,V0] = IntegralAlgo(@matr);
    eigenwerte
    
    for iii=1:length(eigenwerte)
    iii=input('Number:');
    eigenwerte1=eigenwerte(iii);
    V01=V0(:,iii);
    %pause
    yes=1;
    if yes
        reso=100;
        
        x=linspace(-1.6,1.6,reso);
        y=linspace(-1.6,1.6,reso);
        [xx,yy]=meshgrid(x,y);
        z=xx;
        %z=0;
        for i=1:reso
            for j=1:reso
                x=xx(i,j);
                y=yy(i,j);
                [In,~]=inpolygon(x,y,outer(1,:),outer(2,:));
                [In2,~]=inpolygon(x,y,inner(1,:),inner(2,:));
                %In=1;
                if In==1 && In2==0 
                  P=[x;y];
                  summ=0;
                  k=1;
                  L=1;
                  for jj=1:n/2*2    % number of faces
                      % extract the first, middle, and endpoint
                          x1=vi(1,indi1(jj));
                          y1=vi(2,indi1(jj));
                          x2=vi(1,indi2(jj));
                          y2=vi(2,indi2(jj));
                          x3=vi(1,indi3(jj));
                          y3=vi(2,indi3(jj));
                      k=k+2;
                      E1=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),eigenwerte1,1);
                      E2=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),eigenwerte1,2);
                      E3=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),eigenwerte1,3);
              
                      summ=summ+E1*V01(L)+E2*V01(L+1)+E3*V01(L+2);
                      L=L+3;
                  end
                  z(i,j)=summ;
                else
                  z(i,j)=NaN;
                end
            end
        end
       save('complexBump1.mat','z')
       plotthething2('complexBump1.eps',xx,yy,real(z),outer,inner,12,[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2],[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])
    end
    pause
    end
    end
end

function A=matr(k)
  global wavenumber

  wavenumber=k;
  A = createDmatrix();
  A=A+0.5*eye(size(A));
end

function A=createDmatrix()
  global wavenumber
  global vi
  global vhi
  global faces
  global n
  
  indi1=[1:2:n,n+2:2:2*n];
    indi2=[2:2:n,n+3:2:2*n+1];
    indi3=[3:2:n+1,n+4:2:2*(n+1)];
    
  % create the matrix of size n x n
    A=zeros(3*faces*2,3*faces*2);
    M0=zeros(3*faces*2,3*faces*2);
    for i=1:3*faces*2
        % extract collocation point of interval I_i
        xh=vhi(1,i);
        yh=vhi(2,i);
        
        k=1;
        L=1;
        for j=1:n/2*2    % number of faces
            % extract the first, middle, and endpoint
            x1=vi(1,indi1(j));
            y1=vi(2,indi1(j));
            x2=vi(1,indi2(j));
            y2=vi(2,indi2(j));
            x3=vi(1,indi3(j));
            y3=vi(2,indi3(j));
            k=k+2;
            if i==L
              E1=0;
            else
              E1=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,1);
            end
            if i==L+1
              E2=0;
            else
              E2=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,2);
            end
            if i==L+2
              E3=0;
            else
              E3=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,3);
            end
            
            A(i,L  )=E1;
            A(i,L+1)=E2;
            A(i,L+2)=E3;
            if i==L
              E1=0;
            else
              E1=intDL0(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,1);
            end
            if i==L+1
              E2=0;
            else
              E2=intDL0(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,2);
            end
            if i==L+2
              E3=0;
            else
              E3=intDL0(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,3);
            end
            
            M0(i,L  )=E1;
            M0(i,L+1)=E2;
            M0(i,L+2)=E3;
            
            L=L+3;
        end
    end
    
    for i=1:3*faces*2
        A(i,i)=-0.5-sum(M0(i,:),2);
    end
end

function [eigenvalues,Z] = IntegralAlgo(func)
global l
global N
global tol_rank;
global m; 

    ii = 0:1:N;
    tj = 2*pi*ii/N;

    
    while l<=m
        % zufaellige Matrix VDach anlegen
        vD = randn(m,l) +1i*randn(m,l);
        % A0,N und A1,N bestimmen
        A0 = 0;
        A1 = 0;
        for i=1:N
            %i
            p = phi(tj(i));
            ablP = ablPhi(tj(i));
            T=func(p);
            %cond(T)
            A0 = A0 + T\(vD * ablP); % i-ter Summand von A0
            A1 = A1 + T\(vD * p * ablP); % i-ter Summand von A1
        end
        A0 = A0 * 1/(1i*N);
        A1 = A1 * 1/(1i*N);

        % SVD von A0,N
        [V,SW,WH] = svd(A0);
        s = diag(SW);
        s
        k = 0;
        % Rang Test von s
        for i=1:l
            if s(i)> tol_rank
                k=k+1;
            else
                break;
            end
        end
        % falls l=k l erhoehen und wiederholen
        if l==k && ~(l==m)
            l=l+1;
            disp('Ich bin hier')
        else         
            break
        end
    end
    V0 = V(1:m,1:k);
    W0 = WH(1:l,1:k);
    s0 = diag(s(1:k)); 
    % B berechnen mit B = V0
    B = (V0' * A1 * W0) /s0;
    % Eigenwert-Problem fuer B berechnen
    [V,D] = eig(B);
    eigenvalues=diag(D);
    Z=V0*V;
end

% Funktion fuer Phi
function f=phi(t)
global R
global mp
    f =  mp + R * exp(1i*t);
    %f = a * cos(t) + b* sin(t);
end

% Funktion fuer die Ableitung von phi
function f = ablPhi(t)
global R
    f = 1i * R * exp(1i*t);
end

% integrate elastic single layer boundary function
function E=intDL(x1s,x2s,x3s,y1s,y2s,y3s,xhs,yhs,ks,basisf)
  global xh
  global yh
  global x1
  global x2
  global x3
  global y1
  global y2
  global y3
  global basis
  global k
  
  xh=xhs;
  yh=yhs;
  x1=x1s;
  x2=x2s;
  x3=x3s;
  y1=y1s;
  y2=y2s;
  y3=y3s;
  k=ks;
  basis=basisf;

  E=quadgk(@c11_1,0,1);
end

function y=c11_1(s)
  global xh
  global yh
  global x1
  global x2
  global x3
  global y1
  global y2
  global y3
  global k
  global basis
  
    % Jacobian
    dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
    dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
    J=sqrt(dxds.^2+dyds.^2);

    % m_k~
    a=xh-(l1(s)*x1+l2(s)*x2+l3(s)*x3);
    b=yh-(l1(s)*y1+l2(s)*y2+l3(s)*y3);
    r=sqrt(a.^2+b.^2);

    if basis==1
        f=l1a(s);
    elseif basis==2
        f=l2a(s);
    else
        f=l3a(s);
    end

    ny1=dyds./J;
    ny2=-dxds./J;
    
    y=1i*k*besselh(1,k*r)./(4*r).*(a.*ny1+b.*ny2).*J.*f;  % J would cancel out
end
   
% integrate elastic single layer boundary function
function E=intDL0(x1s,x2s,x3s,y1s,y2s,y3s,xhs,yhs,ks,basisf)
  global xh
  global yh
  global x1
  global x2
  global x3
  global y1
  global y2
  global y3
  global basis
  global k
  
  xh=xhs;
  yh=yhs;
  x1=x1s;
  x2=x2s;
  x3=x3s;
  y1=y1s;
  y2=y2s;
  y3=y3s;
  k=ks;
  basis=basisf;

  E=quadgk(@c11_2,0,1);
end

function y=c11_2(s)
  global xh
  global yh
  global x1
  global x2
  global x3
  global y1
  global y2
  global y3
  global basis
  
    % Jacobian
    dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
    dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
    J=sqrt(dxds.^2+dyds.^2);

    % m_k~
    a=xh-(l1(s)*x1+l2(s)*x2+l3(s)*x3);
    b=yh-(l1(s)*y1+l2(s)*y2+l3(s)*y3);
    r=sqrt(a.^2+b.^2);

    if basis==1
        f=l1a(s);
    elseif basis==2
        f=l2a(s);
    else
        f=l3a(s);
    end

    ny1=dyds./J;
    ny2=-dxds./J;
    
    y=1./(2*pi*r.*r).*(a.*ny1+b.*ny2).*J.*f;  % J would cancel out
end
    
    % integrate elastic single layer boundary function
function E=intSL(x1s,x2s,x3s,y1s,y2s,y3s,xhs,yhs,ks,basisf,check)
  global xh
  global yh
  global x1
  global x2
  global x3
  global y1
  global y2
  global y3
  global basis
  global k
  global alpha
  
  xh=xhs;
  yh=yhs;
  x1=x1s;
  x2=x2s;
  x3=x3s;
  y1=y1s;
  y2=y2s;
  y3=y3s;
  k=ks;
  basis=basisf;
 
  if check==0
    E=quadgk(@c11_3,0,1);
  elseif check==1
    E=quadgk(@c11_3,0,alpha)+quadgk(@c11_3,alpha,1);
  elseif check==2
    E=quadgk(@c11_3,0,0.5)+quadgk(@c11_3,0.5,1);
  elseif check==3
    E=quadgk(@c11_3,0,1-alpha)+quadgk(@c11_3,1-alpha,1);
  else
    fprintf('Something is wrong!\n');
  end
end


function y=c11_3(s)
  global xh
  global yh
  global x1
  global x2
  global x3
  global y1
  global y2
  global y3
  global k
  global basis
  
        % Jacobian
        dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
        dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
        J=sqrt(dxds.^2+dyds.^2);
        
        % m_k~
        a=xh-(l1(s)*x1+l2(s)*x2+l3(s)*x3);
        b=yh-(l1(s)*y1+l2(s)*y2+l3(s)*y3);
        r=sqrt(a.^2+b.^2);
        
        if basis==1
            f=l1a(s);
        elseif basis==2
            f=l2a(s);
        else
            f=l3a(s);
        end
        
        y=1i*besselh(0,k*r)/4.*J.*f;
end

% first Lagrange basis function
function y=l1(s)
    u=1-s;
    y=u.*(2*u-1);
end

function y=l1a(s)
    global alpha
    
    u=1-s;
    y=(u-alpha)/(1-2*alpha).*(1-2*s)/(1-2*alpha); % check ok
end

% second Lagrange basis function
function y=l2(s)
    u=1-s;
    y=4*s.*u;
end

function y=l2a(s)
    global alpha
    
    u=1-s;
    y=4*(s-alpha)/(1-2*alpha).*(u-alpha)/(1-2*alpha); % check ok
end

% third Lagrange basis function
function y=l3(s)
    y=s.*(2*s-1);
end

function y=l3a(s)
    global alpha
    
    y=(s-alpha)/(1-2*alpha).*(2*s-1)/(1-2*alpha); % check ok
end


