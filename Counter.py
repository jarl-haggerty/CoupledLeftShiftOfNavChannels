from matplotlib import pyplot
from matplotlib import widgets
import numpy

input = open("experiments/blebbing_1_left_shift_2_stimulation_0.tsv", "r")

names = input.readline().split("\t")
data = []

i = 0

for line in input.readlines():
    data.append(map(float, line.split("\t")))

input.close()

data = numpy.array(data)

print data[1,:]

figure = pyplot.figure()
axis = figure.add_subplot(111)
axis.plot(data[:,0], data[:,1], label = names[1])
axis.plot(data[:,0], data[:,2], label = names[2])
axis.plot(data[:,0], data[:,3], label = names[3])

cursor = widgets.Cursor(axis, useblit = True, color = "black")
cursor.vertOn = False

def onClick(event):
    count_data = data[data[:,0] > data[-1,0]-5000,:]
    print data
    print count_data, len(count_data)

    threshold = cursor.lineh.get_ydata()[0]
    count = 0
    for i in xrange(len(count_data)-1):
        if count_data[i,1] <= threshold and count_data[i+1,1] > threshold:
            count += 1
    print count/5, count

id = figure.canvas.mpl_connect("button_release_event", onClick)
print numpy.nonzero(data[:,0] < data[-1,0]-5000)
pyplot.show()
