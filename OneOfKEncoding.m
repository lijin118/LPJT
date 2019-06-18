function Yk = OneOfKEncoding(Y)
Yk = full(sparse(1:length(Y), Y, ones(length(Y),1)));
end
