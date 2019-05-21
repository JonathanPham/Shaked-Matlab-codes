import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import string
import sys

#input parameter check
if len(sys.argv) < 3:
    print('Usage:')
    print('python3 postprocess.py <inputfile> <solutionfile>')
    sys.exit(0)

#reading data
inpf=open(sys.argv[1],'r')
solf=sys.argv[2]
lines=inpf.readlines()
length=float(lines[0].split()[0])
width=float(lines[0].split()[1])
h=float(lines[0].split()[2])
low=float(lines[1].split()[0])
high=float(lines[1].split()[1])
cols=int(length/h)
xlist=[]
count=0
with open(solf) as f:
    for line in f:
        if (count%cols==0):
            left=float(line)
        xlist.append(float(line))
        if (count%cols==cols-1):
            xlist.append(left)
        count+=1
meant=np.mean(xlist)
#reshaping data
xlist.reverse()
temparray=np.array(xlist)
x=np.arange(0,length+h,h)
y=np.arange(0,width+h,h)
temp2D=temparray.reshape(len(y),len(x))
# print(temp2D[0,:])
X, Y =np.meshgrid(x,y)

#isoline
isoline=np.zeros(temp2D.shape[1])
for ind in range(temp2D.shape[1]):
    isoline[ind]=np.interp(meant,temp2D[:,ind],Y[:,ind])

#Console output
print("Input file processed: ",  sys.argv[1])
print("Mean Temperature: ",  meant)

#plotting everything
fig = plt.figure()
plt.pcolor(X, Y, temp2D, cmap='jet')
plt.plot(x, isoline)
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Temperature Distribution in the Pipe')
plt.axis([0, length, -1*width, 2*width])
plt.savefig(solf + ".png")

