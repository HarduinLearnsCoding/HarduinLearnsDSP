function [out] = CONJ32(in)
re=real(in);
img=-imag(in);
out=complex(re,img);
out=int32(out);
end

