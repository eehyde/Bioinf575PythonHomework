import sys

fileName = sys.argv[1]
openedfile = open(fileName,"r")

numlines = 0
numwords = 0
numchars = 0
for x in openedfile:
	words = x.split()
	numlines += 1
	numwords += len(words)
	numchars += len(x)

print(numlines, numwords, numchars)	






