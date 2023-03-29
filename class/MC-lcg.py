x0=int(input("seed="))

x=(5*x0+3)%8; print(x)
while x != x0 :
	x=(5*x+3)%8
	print(x)
