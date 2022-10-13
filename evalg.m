function [g] = evalg(X)

global A B

g = A' * ( A * X - B );