function out = steer(n, d, phi);
out = exp(-1j*2*pi*([0:n-1].')*d*sind(phi));
end

