function [out] = CONJ16(in)
re=real(in);
img=-imag(in);
out=complex(re,img);
out=int16(out);
end

