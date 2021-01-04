function plotit()
  f=@(R)besseljp(0,R);
  fzero(f,4)
  f=@(R)besseljp(1,R);
  fzero(f,3)
  
  g=@(R)besselj(1,1.841184*R);
  
  reso=100;
        
  x=linspace(-2.1,2.1,reso);
  y=linspace(-2.1,2.1,reso);
  [xx,yy]=meshgrid(x,y);
  z=xx;
  for i=1:reso
    for j=1:reso
      x=xx(i,j);
      y=yy(i,j);
      [theta,R]=cart2pol(x,y);
      z(i,j)=g(R)*sin(theta);
    end
  end
  
  n=100;
  ii=0:1:n;
    ti=2*pi*ii/n;
    vi=[cos(ti);sin(ti)];
    outer=vi(:,1:n+1);
    
  figure(1)
  contourf(xx,yy,z,30);
  line(outer(1,:),outer(2,:),'Color','k','Linewidth',2)
  colorbar()
end

function z=besseljp(n,x)
  z=n./x.*besselj(n,x)-besselj(n+1,x);
end