# -*- octave -*-
function td1_exo2 (i)

  ## usage:  td1_exo2 (N)
  ##
  ##
  eval(sprintf('temperature%03d', i));
  T=eval(sprintf('T%d', i));
  N=size(T,1);
  [xx,yy]=meshgrid(linspace(0,1,N),linspace(0,1,N));
  contour(xx,yy,T);
  #mesh(xx,yy,T); 
endfunction
