function [result] = SUB16_CX(input1,input2)
re=int32(int16(real(input1))-int16(real(input2)));
im=int32(int16(imag(input1))-int16(imag(input2)));
result=complex(re,im);
end
