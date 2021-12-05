function [result] = ADD32_CX(input1,input2)
re=int32(double(real(input1))+double(real(input2)));
img=int32(double(imag(input1))+double(imag(input2)));
result=complex(re,img);
end

