function [h]=flat_pcolor(x, y, c);
% FLAT_PCOLOR makes a pseudocolor plot of the entire matrix c
%
%call flat_pcolor(c) or flat_pcolor(x, y, c)
%
%FLAT_PCOLOR makes a pseudocolor plot of the entire matrix c using the
%same input arguments as pcolor but does so in the way that one would
%intuitively expect pcolor to work.  The problem with pcolor is that if
%shading is 'flat' or 'faceted' then the value of c(i,j) is used to draw
%a rectangle where the bottom left hand corner is at (x(j), y(i)) and
%the top right hand corner is at (x(j+1), y(i+1)).  Accordingly, the
%last row of c and last column of c are not used.  flat_pcolor uses x
%and y (which may be vectors or matrices) as the center of the rectangle
%and all rows and columns of c are used.  If only one argument is passed
%to flat_pcolor then this is assumed to be the matrix c and x and y are
%assumed to be the column and row indices.  The shading is set to
%'flat'.  Note that if shading is 'interp' then pcolor does a slow,
%bilinear interpolation between data points as would be expected but
%this means that a single NaN creates 4 blank spaces.  Since flat_pcolor
%shifts the x and y axes by half a grid length then you should not set
%shading('interp') after running flat_pcolor - use pcolor instead.  Note
%also that if x or y is a scalar then an axis going from 0 to 1 will be
%used.  This is useful for drawing a legend.
%
% CALLER:   general purpose
% CALLEE:   pcolor
% VERSION 1.1
%     Copyright J. V. Mansbridge, CSIRO, june 28 1993

% Create interpolated and extrapolated vectors or matrices in the x
% and y directions.

if nargin == 1;
  
  c = x;
  [m, n] = size(c);
  xs = 0.5:(n + 0.5);
  ys = 0.5:(m + 0.5);
  
elseif nargin == 3;
  
  [m, n] = size(c);

  % First the x vector or matrix.

  [mx, nx] = size(x);
  if nx == 1;
    if mx == 1;    % scalar case
      xs = [ 0 1 ];
    else          % column vector case
      x0 = 2*x(1) - x(2);
      xf = 2*x(mx) - x(mx-1);
      xs = ([ x0; x ] + [ x; xf ])/2;
    end
  else
    if mx == 1;   % row vector case
      x0 = 2*x(1) - x(2);
      xf = 2*x(nx) - x(nx-1);
      xs = ([ x0 x ] + [ x xf ])/2;
    else          % matrix case
      x0 = 2*x(1, :) - x(2, :);
      xf = 2*x(mx, :) - x((mx-1), :);
      xs = ([ x0; x ] + [ x; xf ])/2;
      x0 = 2*xs(:, 1) - xs(:, 2);
      xf = 2*xs(:, nx) - xs(:, (nx-1));
      xs = ([ x0 xs ] + [ xs xf ])/2;
    end
    
  end

  % Now the y vector or matrix.
  
  [my, ny] = size(y);
  if ny == 1;
    if my == 1;    % scalar case
      ys = [ 0 1 ];
    else          % column vector case
      y0 = 2*y(1) - y(2);
      yf = 2*y(my) - y(my-1);
      ys = ([ y0; y ] + [ y; yf ])/2;
    end
  else
    if my == 1;    % row vector case
      y0 = 2*y(1) - y(2);
      yf = 2*y(ny) - y(ny-1);
      ys = ([ y0 y ] + [ y yf ])/2;
    else          % matrix case
      y0 = 2*y(1, :) - y(2, :);
      yf = 2*y(my, :) - y((my-1), :);
      ys = ([ y0; y ] + [ y; yf ])/2;
      y0 = 2*ys(:, 1) - ys(:, 2);
      yf = 2*ys(:, ny) - ys(:, (ny-1));
      ys = ([ y0 ys ] + [ ys yf ])/2;
    end
  end

end

% Create an extended version of the c matrix with an extra row and
% column at the ends.

cs    = [ c      , c(:,n);
          c(m,:) , c(m, n)  ];

% Call pcolor to produce the picture.

[h]=pcolor(xs, ys, cs);
shading('flat');

