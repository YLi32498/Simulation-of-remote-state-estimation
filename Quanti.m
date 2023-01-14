function S_q = Quanti(S,C,L,N,Index)
if Index == 0; % encoding
    for i = 1:length(S)
        S_q(1,i) = floor(  (S(i)-C(i)+L(i))*N./L(i)/2  );
    end
else if Index == 1 % decoding
        for i = 1:length(S)
            S_q(1,i) = S(i)*2*L(i)./N + C(i) -L(i) + L(i)./N;
        end
% else if Index == 2 % encoding + decoding
%         for i = 1:length(S)
%             S_e(1,i) = floor(  (S(i)-C(i)+L(i)./2)*N./L(i)  );
%             S_q(1,i) = S_e(i)*L(i)./N + C(i) -L(i)/2;
%         end
else 
    error('wrong Index number')
end
end