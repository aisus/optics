function r=swap_parts(arr)
  n = length(arr);
  half = floor(n / 2);
  r = vertcat(arr(half + 1:n), arr(1:half));
endfunction
