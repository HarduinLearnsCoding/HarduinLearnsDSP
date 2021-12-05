function [output_vector] = Interleave(fixed_input)

output_vector=zeros(2,length(fixed_input));

%Need to fix according to dimensions of fixed_input

for i=1:length(fixed_input)
    output_vector(1,i)=int16(real(fixed_input(i)));
    output_vector(2,i)=int16(imag(fixed_input(i)));
end

end

