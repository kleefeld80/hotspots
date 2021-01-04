function err=HelmDL2(n)
    % parameters
    wavenumber=2;
    
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
    
    % calculate normal in the collocation point exactly
    % r(t) -> T(t)=r'(t)/|r'(t)| -> N(t)=T'(t)/|T'(t)|
    % Sphere
    % vec=[cos(si(i));sin(si(i))];
    % Ellipsoid
    % vec=[a*b^2*cos(si(i));a^2*b*sin(si(i))];
    %nx=zeros(2,n);
    for i=1:n
        vec=[a*b^2*cos(si(i));a^2*b*sin(si(i))]; % Circle and Ellipse
        nx(:,i)=vec/norm(vec,2);
    end
    
    % calculate approximate normal for each segment
    %ny=zeros(2,n);
    %k=1;
%     for i=1:n/2
%         % calculate slope of line segment
%         indi1=k;
%         indi2=piminus(k,n);
%         indi3=piminus(k+1,n);
%         x1=vi(1,indi1);
%         y1=vi(2,indi1);
%         x2=vi(1,indi2);
%         y2=vi(2,indi2);
%         x3=vi(1,indi3);
%         y3=vi(2,indi3);
%         
%         s=0;
%         dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
%         dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
%         J=sqrt(dxds.^2+dyds.^2);
%     
%         ny(:,k)  =[dyds;-dxds]/J;
%         
%         s=0.5;
%         dxds=(4*s-3)*x1+4*(1-2*s)*x2+(4*s-1)*x3;
%         dyds=(4*s-3)*y1+4*(1-2*s)*y2+(4*s-1)*y3;
%         J=sqrt(dxds.^2+dyds.^2);
%         ny(:,k+1)=[dyds;-dxds]/J;
%         k=k+2;
%     end
    ny=nx;
    hold on
    for i=1:n
        plot([0,ny(1,i)]+vhi(1,i),[0,ny(2,i)]+vhi(2,i))
    end
    hold off
    %nx
    %ny
    %nx-ny
    
    % create the right hand size (boundary condition)
    rhs=zeros(n,1);
    for i=1:n
        % extract collocation point of interval I_i
        xh=vhi(1,i);
        yh=vhi(2,i);
        rhs(i,1)=truesolution(xh,yh,wavenumber);
    end
    
    % create the matrix of size n x n
    A=zeros(n,n);
    for i=1:n
        % extract collocation point of interval I_i
        xh=vhi(1,i);
        yh=vhi(2,i);
        
        k=1;
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
            E1=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,1);
            E2=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,2);
            E3=intDL(x1,x2,x3,y1,y2,y3,xh,yh,wavenumber,3);
            A(i,indi1)=A(i,indi1)+E1;
            A(i,indi2)=A(i,indi2)+E2;
            A(i,indi3)=A(i,indi3)+E3;
        end
    end
    
    %A
    
    lsg=(0.5*eye(n)+A)\rhs;
    
    % calculate elastic solution in the exterior
    P=[5;5];
    summ=0;
    k=1;
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

        summ=summ+E1*lsg(indi1)+E2*lsg(indi2)+E3*lsg(indi3);
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

  E=quadgk(@c11_1,0,1,'AbsTol',1e-7);
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
        f=l1(s);
    elseif basis==2
        f=l2(s);
    else
        f=l3(s);
    end

    ny1=dyds./J;
    ny2=-dxds./J;
    
    y=1i*k*besselh(1,k*r)./(4*r).*(a.*ny1+b.*ny2).*J.*f;  % J would cancel out
end
    
% first Lagrange basis function
function y=l1(s)
    u=1-s;
    y=u.*(2*u-1);
end

% second Lagrange basis function
function y=l2(s)
    u=1-s;
    y=4*s.*u;
end

% third Lagrange basis function
function y=l3(s)
    y=s.*(2*s-1);
end

% true solution, first column of fundamental matrix
function f=truesolution(a,b,k)
      r=sqrt(a.^2+b.^2);
      f=1i*besselh(0,k*r)/4;
end

