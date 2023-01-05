function [ACorrected1,CCorrected1,P]=PressureSolve(A,C)

[Asub,idx]=licols(A,1e-10);
[size1,size2]=size(idx);
for iii=1:size2
    ACorrected1(iii,:)=Asub(idx(iii),:);
    CCorrected1(iii,1)=C(idx(iii),1);
end

PCorrected1=inv(ACorrected1)*CCorrected1;

P=zeros(size(C));
for iii=1:size2
    P(idx(iii),1)=PCorrected1(iii);
end

end
