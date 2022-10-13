function [f] = evalf(X)

global A B

f = 0.5 * norm( A * X - B,'fro')^2;