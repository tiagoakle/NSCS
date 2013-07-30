function [F,G,H,feas] = BarrFuncPowHild(x1,x2,x3,K,want)

for j = 1:size(x1,1)
	[FF,FF1,FF2,Ffeas] = EHbarrier(x2(j),x3(j),x1(j),K.powc.alph(j))
end