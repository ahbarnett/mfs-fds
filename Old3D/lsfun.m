

function x = lsfun(b)
global r A
x = r\(r'\(A'*b));
x = x + r\(r'\(A'*(b - A*x)));
end
