function [Ne dNe]=buildnlin_bound(s,r)
Ne=1/4*[(1-s)*(1-r), (1+s)*(1-r), (1+s)*(1+r), (1-s)*(1+r)];
dNe=1/4*[-(1-r), (1-r), (1+r), -(1+r); -(1-s), -(1+s), (1+s), (1-s)];
end