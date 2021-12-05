function [result] = SUB32_CX(input1,input2)
re=int64(int32(real(input1))-int32(real(input2)));
img=int64(int32(imag(input1))-int32(imag(input2)));
result=complex(re,img);
end

