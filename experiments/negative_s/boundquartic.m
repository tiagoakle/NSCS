samples = 1000
x = [0:samples-1]./(samples*100)

f = @(l) -l^2/(1-l) + 1/(1-l)

y = arrayfun(f,x)

plot(y,x)