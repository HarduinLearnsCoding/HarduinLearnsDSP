function [fixed_output,max_norm] = ADConvert(float_input,max_norm)

if nargin==1
    max_norm=max((abs(float_input).^2),[],'all');
end

%Modify based on the characteristics of the frame float_input

A=8-1;

fixed_output=(A*float_input/max_norm);

fixed_output=int16(fixed_output);

end

