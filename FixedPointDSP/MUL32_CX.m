function [result] = MUL32_CX(input1,input2)
r1=real(input1);
r2=real(input2);
i1=imag(input1);
i2=imag(input2);

re=int64(int32(double(r1)*double(r2))-int32(double(i1)*double(i2)));
img=int64(int32(double(r1)*double(i2))+int32(double(i1)*double(r2)));
result=complex(re,img);
end

