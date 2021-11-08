clc;
size=128;
dftnew=fft(eye(size));
[createdarray]=ffttable(128,1,128,1,128);
% % Checking the equivalence
% count=0;
% for i=1:size
%     for j=1:size
%         if createdarray~=dftnew
%             count=count+1;
%         end
%     end
% end
% if count==0
%     disp('yes');
% end
