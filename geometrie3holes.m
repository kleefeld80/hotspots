function points=geometrie3holes(n)
  eps=1/2;
  
  % Part 1
  points1=ellipsepart(4,8*eps,1,2*eps,-pi/2,-pi,n,1);
  points1=[points1 ellipsepart(2,8*eps,1,2*eps,0,pi,n,1)];
  points1=[points1 ellipsepart(0,8*eps,1,eps,0,-pi,n,1)];
  points1=[points1 ellipsepart(-2,8*eps,1,2*eps,0,pi,n,1)];
  points1=[points1 ellipsepart(-4,8*eps,1,2*eps,0,-pi/2,n,0)];
  %pointsrand=ellipsepart(-4,0,8*eps,6*eps,pi/2,3*pi/2,n,1);     % 8->6
  points2=fliplr(points1);
  points1=ellipsepart(4,8*eps,1,2*eps,-pi/2,-pi,n,1);
  points1=[points1 ellipsepart(2,8*eps,1,2*eps,0,pi,n,1)];
  points1=[points1 ellipsepart(0,8*eps,1,eps,0,-pi,n,1)];
  points1=[points1 ellipsepart(-2,8*eps,1,2*eps,0,pi,n,1)];
  points1=[points1 ellipsepart(-4,8*eps,1,2*eps,0,-pi/2,n,0)];
  pointsrand=ellipsepart(-4,0,6*eps,6*eps,pi/2,3*pi/2,n,1);          % 3. Arg 6->8
  points2(2,:)=-points2(2,:);
  pointsrand2=ellipsepart(4,0,6*eps,6*eps,-pi/2,pi/2,n,0);
  pointsp1=[points1(:,1:end-1) pointsrand points2(:,1:end-1) pointsrand2];
  
  % Part2
  points1=ellipsepart(4,3*eps,1,2*eps,pi/2,pi,n,1);
  points1=[points1 ellipsepart(2,3*eps,1,2*eps,0,-pi,n,1)];
  points1=[points1 ellipsepart(0,3*eps,1,eps,0,pi,n,1)];
  points1=[points1 ellipsepart(-2,3*eps,1,2*eps,0,-pi,n,1)];
  points1=[points1 ellipsepart(-4,3*eps,1,2*eps,0,pi/2,n,0)];
  pointsrand=ellipsepart(-4,0,5*eps,5*eps,pi/2,3*pi/2,n,1);     % 7->5
  points2=fliplr(points1);
  points2(2,:)=-points2(2,:);
  pointsrand2=ellipsepart(4,0,5*eps,5*eps,-pi/2,pi/2,n,0);
  pointsp2=fliplr([points1(:,1:end-1) pointsrand points2(:,1:end-1) pointsrand2]);
  
  % Circle
  m=length(pointsp1)-1;
  ii=0:1:m;
  ti=2*pi*ii/m;
  R=0.8
  mux=2;
  muy=2.75;
  pointsp3=fliplr([mux+R*cos(ti);muy+R*sin(ti)]);
  
  % Circle
  m=length(pointsp1)-1;
  ii=0:1:m;
  ti=2*pi*ii/m;
  mux=-2;
  muy=2.75;
  pointsp4=fliplr([mux+R*cos(ti);muy+R*sin(ti)]);
  
  points=[pointsp1 pointsp2 pointsp3 pointsp4];
  
  %points=[pointsp1];
  %figure(1)
  %clf
  %hold on
  %fill(pointsp1(1,:),pointsp1(2,:),[0.5 0.5 0.5])
  %fill(pointsp2(1,:),pointsp2(2,:),'w')
  %axis equal
  %hold off
end

% creates part of an ellipse
function points=ellipsepart(x0,y0,a,b,phi0,phi1,n,erase)
  phi=linspace(phi0,phi1,n);
  x=x0+a*cos(phi);
  y=y0+b*sin(phi);
  
  points=[x;y];
  if erase
    points=points(:,1:end-1);
  end
end