%% -*- Mode: octave -*-
function polyplot (fname)

  ## usage:  polyplot (fname)
  ##
hold off;
M=load(fname);
for i=2:size(M,2)
  plot(M(:,1),M(:,i),[';P_{' int2str(i-2) '};']);
  hold on
end
hold off
endfunction


