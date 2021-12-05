function [result] = SHIFT_LONG(input1,input2)
result=int64(bitsra(int64(input1),int16(input2)));
end

