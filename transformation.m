function [ MB ] = transformation( MA, vv, tv,insv)
%Summary of this function goes here
% function that transforms the matrix A according to the 
% transformation vv, tv and iv
%   Detailed explanation goes here

%opeartions
% ih iv dh dv
ih=1;
iv=2;
dh=3;
dv=4;


%sort the vv to begin
[vv(:,1) index ]=sort(vv(:,1));
vv(:,2)=vv(index,2);

%initialize the matrix result matrix
MB=MA;

%keep a count of the number of insertions
Nins=1;


for  i1=1:length(vv)        
            
        switch  tv(i1)
        
        case ih 
               %shift the row out
               MB(vv(i1,1), vv(i1,2)+1: end+1)=MB(vv(i1,1), vv(i1,2): end);
               %insert the element
               MB(vv(i1))=insv(Nins); 
               Nins=Nins+1;
               
               %update tranformation vector
               %shift transformations on the same row by one
               index=find(vv(:,2)>vv(i1,1));
               vv(index,2)=vv(index,2)+1;
               
        case iv 
               %shift the column out
               MB(vv(i1,2)+1: end+1, vv(i1,2))=MB(vv(i1,2): end, vv(i1,2));
               %insert the element
               MB(vv(i1))=insv(Nins); 
               Nins=Nins+1;
               
               %update tranformation vector
               %shift transformations on the same row by one
               index=find(vv(:,1)>vv(i1,1));
               vv(index,1)=vv(index,1)+1;
               
               
        case dh
               %shift the row in
               MB(vv(i1,1), vv(i1,2)-1: end-1)=MB(vv(i1,1), vv(i1,2): end);               

               %update tranformation vector
               %shift transformations on the same row by one
               index=find(vv(:,1)>=vv(i1,1));
               vv(index,1)=vv(index,1)-1;
        case dv
               %shift the column in
               MB(vv(i1,2)-1: end-1, vv(i1,2))=MB(vv(i1,2): end, vv(i1,2));                                                
               
               %update tranformation vector
               %shift transformations on the same row by one
               index=find(vv(:,2)>=vv(i1,1));
               vv(index,2)=vv(index,2)+1;
        end

end

%check if you can remove zeros.

MB=MB(:,sum(MB)~=0);

MB=MB(sum(MB')~=0,:);

%NOTE: is there a way of reducing the complexity of this?
% the commented code is about that. 

% 
% %initialize the 
% currentPos=[0 0];
% 
% 
% for  i1=1:length(vv)
%     
%     % previous transformation 
%     prevPov=currentPos
%     
%     % current transformation
%     currentPos=vv(i1)
%         
%     %start by copying the matrix between two successive transformations
%     [r c]=find(BM==0);    
%             
%     % copy the rows
%     
%     %find the edges of the current matrix
%     for i2=1:prevPov(1) 
%         cE(i2)=max(c(find(r=i2)))        
%     end
%         
%     columLength=prevPov(2)-currentPos(2)
%     
%     for i2=1:prevPov(1)
%         MB(i2;cE(i2):cE(i2)+columLength)=MA(i2; prevPov(2):prevPov(2)+columLength)
%     end
%         
%     % copy the columns
%     
%     %find the edges of the current matrix
%     for i2=1:prevPov(2) 
%         vE(i2)=max(v(find(c=i2)))        
%     end
%         
%     rowLength=prevPov(1)-currentPos(1)
%     
%     for i2=1:prevPov(1)
%         MB(i2;cE(i2):cE(i2)+columLength)=MA(i2; prevPov(2):prevPov(2)+columLength)
%     end
%         
%     
%     % copy the remaining ones
%     
%     
%     MB() =MA(prevPov(1)-currentPos(1), prevPov(2)-currentPos(2))
%     
%     
%     typeT=tv(i1);
%     
%     switch  typeT
%         
%         case ih 
%             
%             
%             
%         case iv 
%         
%         case dh
%             
%         case dv
%            
%             
% 
%         
% end