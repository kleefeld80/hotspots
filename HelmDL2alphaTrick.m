% tested successful
% alpha=0.1 and 0.2 worse convergence rates
% alpha=(1-sqrt(3/5))/2 rate=3
function err=HelmDL2alphaTrick(n)
    global alpha
    
    alpha=(1-sqrt(3/5))/2;
    faces=n/2;
    fprintf('Number of faces: %d\n',faces);
    fprintf('Number of collocation nodes: %d\n',3*faces);
    
    % parameters
    wavenumber=1;
    
    a=1.2;
    b=1;
    
    % create triangulation mit knodes v_i
    ii=0:1:n;
    ti=2*pi*ii/n;
    vi=[a*cos(ti);b*sin(ti)]; % Circle and Ellipse
    %vi=[1.5*sin(ti);cos(ti)+0.65*cos(2*ti)-0.65]; % Kite
    %vi=[sqrt(cos(ti).*cos(ti)+4*sin(ti).*sin(ti)).*cos(ti);sqrt(cos(ti).*cos(ti)+4*sin(ti).*sin(ti)).*sin(ti)]; % Peanut
    
    % create collocation nodes v^_i
    vhi=vi;
    si=ti;
    %for i=1:n
    %    x1=vi(1,i);
    %    y1=vi(2,i);
    %    x2=vi(1,i+1);
    %    y2=vi(2,i+1);
    %    vhi(:,i)=[x1+x2;y1+y2]/2;
    %end
    %vhi=[1.5*sin(si);cos(si)+0.65*cos(2*si)-0.65];
    %vhi=[sqrt(cos(si).*cos(si)+4*sin(si).*sin(si)).*cos(si);sqrt(cos(si).*cos(si)+4*sin(si).*sin(si)).*sin(si)]; % Peanut
    
    hold on
    plot(vi(1,:),vi(2,:));
    plot(vhi(1,:),vhi(2,:),'x');
    hold off
    
    vhi=vi;
    k=1;
    L=1;
    for i=1:faces
      indi1=k;
      indi2=piminus(k,n);
      indi3=piminus(k+1,n);
      x1=vi(1,indi1);
      y1=vi(2,indi1);
      x2=vi(1,indi2);
      y2=vi(2,indi2);
      x3=vi(1,indi3);
      y3=vi(2,indi3);
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
    
    % create the right hand size (boundary condition)
    rhs=zeros(3*faces,1);
    for i=1:3*faces
        % extract collocation point of interval I_i
        xh=vhi(1,i);
        yh=vhi(2,i);
        rhs(i,1)=truesolution(xh,yh,wavenumber);
    end
    
    % create the matrix of size n x n
    A=zeros(3*faces,3*faces);
    M0=zeros(3*faces,3*faces);
    for i=1:3*faces
        % extract collocation point of interval I_i
        xh=vhi(1,i);
        yh=vhi(2,i);
        
        k=1;
        L=1;
        for j=1:n/2    % number of faces
            % extract the first, middle, and endpoint
            indi1=k;
            indi2=piminus(k,n);
            indi3=piminus(k+1,n);
            x1=vi(1,indi1);
            y1=vi(2,indi1);
            x2=vi(1,indi2);
            y2=vi(2,indi2);
            x3=vi(1,indi3);
            y3=vi(2,indi3);
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
    
    for i=1:3*faces
        A(i,i)=-0.5-sum(M0(i,:),2);
    end
    
    lsg=(0.5*eye(3*faces)+A)\rhs;
    
    % calculate elastic solution in the exterior
    P=[5;5];
    summ=0;
    k=1;
    L=1;
    for j=1:n/2    % number of faces
        % extract the first, middle, and endpoint
        indi1=k;
        indi2=piminus(k,n);
        indi3=piminus(k+1,n);
        x1=vi(1,indi1);
        y1=vi(2,indi1);
        x2=vi(1,indi2);
        y2=vi(2,indi2);
        x3=vi(1,indi3);
        y3=vi(2,indi3);
        k=k+2;
        E1=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),wavenumber,1);
        E2=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),wavenumber,2);
        E3=intDL(x1,x2,x3,y1,y2,y3,P(1),P(2),wavenumber,3);

        summ=summ+E1*lsg(L)+E2*lsg(L+1)+E3*lsg(L+2);
        L=L+3;
    end
    
    sol=summ;
    tru=truesolution(P(1),P(2),wavenumber);
    err=norm(abs(sol-tru),2);
    fprintf('Approx. %f %f\n',real(sol),imag(sol));
    fprintf('True    %f %f\n',real(tru),imag(tru));
    fprintf('Error.  %g\n',err);
end

function indi=piminus(j,n)
    if j==n
        indi=1;
    else
        indi=j+1;
    end
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

% true solution, first column of fundamental matrix
function f=truesolution(a,b,k)
      r=sqrt(a.^2+b.^2);
      f=1i*besselh(0,k*r)/4;
end

