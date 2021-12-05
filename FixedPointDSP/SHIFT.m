function [result] = SHIFT(input1,input2)
result=int32(bitsra(int32(input1),int16(input2)));
end

