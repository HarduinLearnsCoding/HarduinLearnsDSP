function [result] = MUL16_CX_old(input1,input2)
r1=real(input1);
r2=real(input2);
i1=imag(input1);
i2=imag(input2);

re=int32(int32(int16(r1)*int16(r2))-int32(int16(i1)*int16(i2)));
img=int32(int32(int16(r1)*int16(i2))+int32(int16(i1)*int16(r2)));
result=complex(re,img);
end

