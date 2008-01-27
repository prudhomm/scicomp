%% -*- Mode: octave -*-

hold off;
M=load('legendre.dat');
for i=2:size(M,2)
  plot(M(:,1),M(:,i),[';L_{' int2str(i-2) '};']);
  hold on
end
hold off