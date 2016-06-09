function out = testsimp2D(table,a,b,c,d,NX,NY)
% clc; clear;
% a = 0; b = 1;
% c = a; d = b;
% NX = 512; NY = NX;
% x = linspace(a,b,NX+1);
% y = x;
% func = @(x,y) cos(x.^2 + y.^2);
% table = func(x,y);
% for i = 1:length(x)
% for j = 1:length(y)
% table(i,j) = func(x(i),y(j));
% end
% end
% break;


% Ensure the number of intervals is even!
NX = 2*ceil(NX/2);
NY = 2*ceil(NY/2);

% Set up the integration step sizes in the x and y directions
hx = (b - a)/NX;
hy = (d - c)/NY;
% break;
% define grid vectors
% yg = c:hy:d;
% xg = a:hx:b;
% [xxg, yyg] = ndgrid(xg,yg);

% Note that the ndgrid routine is only included inside simp2D to make it
% self-contained. If you were using simp inside a loop, it would be better to
% calculate the vectors and ndgrid step outside the function and pass xxg, yyg
% and a, b, c, d as globals say. This can save 25% in run time in large loops.
% For single use or small loops the penalty is small.

% Now set up a matrix U that contains the values of the function evaluated at all
% points on the 2D grid setup by xg and yg.

% U = feval(func,xxg,yyg);
U = table;

% Evaluate the contribution from the corner points first.
% These all have weight 1. NB U(1,1) corresponds to func(a,b) etc.

s1 = ( U(1,1) + U(1,NY+1) + U(NX+1,1) + U(NX+1,NY+1) );

% Now sum the contributions from the terms along each edge not including
% corners. There are 4 edges in the 2D case that contribute to the sum
% and we have points with weight 4 and points with weight 2. Points
% with weight 4 are acessed by indices 2:2:N (N=NX,NY,NZ), while points with
% weight 2 are accessed by indices 3:2:N-1.

% Define vectors of odd and even indices for each direction:

ixo = 2:2:NX;
ixe = 3:2:NX-1;
iyo = 2:2:NY;
iye = 3:2:NY-1;

s2 = 2*( sum(U(1,iye)) + sum(U(NX+1,iye)) + sum(U(ixe,1)) + sum(U(ixe,NY+1)) );
s3 = 4*( sum(U(1,iyo)) + sum(U(NX+1,iyo)) + sum(U(ixo,1)) + sum(U(ixo,NY+1)) );


% Now we look at the remaining contributions on the interior grid points. 
% Looking at our array example above we see that there
% are only 3 different weights viz. 16, 8 and 4. Some thought will show that
% using our definitions above for odd and even gridpoints, that weight 16 is
% only found at points (xodd, yodd), weight 4 is found at points (xeven,yeven)
% while weight 8 is found at both (xodd,yeven) or (xeven,yodd).


% Our contribution from interior points is then

s4 = 16*sum( sum( U(ixo,iyo) ) ) + 4*sum( sum( U(ixe,iye) ) );
s5 =  8*sum( sum( U(ixe,iyo) ) ) + 8*sum( sum( U(ixo,iye) ) );

% Finally add all the contributions and multiply by the step sizes hx, hy and 
% a factor 1/9 (1/3 in each direction).


out = s1 + s2 + s3 + s4 + s5;
out = out*hx*hy/9.0;