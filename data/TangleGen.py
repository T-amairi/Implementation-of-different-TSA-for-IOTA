#modules
from graphviz import Digraph 
from PIL import Image
import os
import csv

#paths
infile = r'C:\Users\AMAIRI\Desktop\MAIN\IOTA\data\Tracking\TrackerTangle.txt'
outfile = r'C:\Users\AMAIRI\Desktop\MAIN\IOTA\data\Tangle.txt'
path = r'C:\Users\AMAIRI\Desktop\MAIN\IOTA\data'

#filter 
if os.path.exists(infile):
    s = set()
    with open(outfile, 'w') as out:
        for line in open(infile):
            if not line.strip(): continue
            if line not in s:
                out.write(line)
                s.add(line)
    os.remove(infile)
    out.close()

#dict for sites with key = time of sim and val = sites
sites = dict()
temp = 0

#create the dict 
with open(outfile, 'r',) as file:
    reader = csv.reader(file, delimiter = ';')
    for row in reader:
        time = row[0]
        site_name = row[1]
        site_neib = row[2]
        if temp != time:
            temp = row[0]
            sites[time] = []
            sites[time].append((site_name,site_neib.split(",")))
        else:
            sites[time].append((site_name,site_neib.split(",")))

#compute the number of tips and return the list of Tips 
def Tips(listSite):
    numberTips = 0
    listTips = []
    for node in listSite:
        if node[1][0] == '':
            listTips.append(node[0])
            numberTips += 1
    return numberTips,listTips

#compute the number of tips that approve each sites 
def TipsApp(listSite,listTips):
    listSite.reverse()
    dictApp = {}
    for node in listSite:
        dictApp[node[0]] = []
    for node in listSite:
        if not node[0] in listTips : 
            for neib in node[1]:
                if neib in listTips:
                    dictApp[node[0]].append(neib)
                else:    
                    dictApp[node[0]] = list(set(dictApp[node[0]] + dictApp[neib]))
    return dictApp
        
#to get maxTime of the sim 
maxTime = 0.0
for t, s in sites.items():
    if maxTime < float(t):
        maxTime = float(t)

#colors of the node for graphviz 
colors = ['green','orange','red']
#render each Tangle of each time sim if 1, the last Tangle otherwise if 0 
test = 0
#confidence 
conf = 0.5

#generating the Tangle 
if test:
    for t, s in sites.items():
        g = Digraph(comment='Tangle' + t + 'sec',format='png',node_attr={'shape': 'box','style': 'filled'})
        g.attr(rankdir='LR')
        g.attr(label=r'Tangle at ' + t + ' sec (simulation time)',labelloc='t')

        numberTips,listTips = Tips(s)
        dictApp = TipsApp(s,listTips)

        for node in s:
            if len(dictApp[node[0]]) >= int(numberTips*conf):
                g.node(node[0],color=colors[0])
            elif node[1][0] != '':
                g.node(node[0],color=colors[1])
            else:
                g.node(node[0],color=colors[2])

        for node in s:
                for neib in node[1]:
                    if neib != '':
                        g.edge(node[0],neib,dir="back")
        g.render(path +'\im\Tangle' + t + 'sec')
        os.remove(path +'\im\Tangle' + t + 'sec')

else:
    t = str(maxTime)
    s = sites[t]
    numberTips,listTips = Tips(s)
    dictApp = TipsApp(s,listTips)
    g = Digraph(comment='Tangle' + t + 'sec',format='png',node_attr={'shape': 'box','style': 'filled'})
    g.attr(rankdir='LR')
    g.attr(label=r'Tangle at ' + t + ' sec (simulation time)',labelloc='t')

    for node in s:
        if len(dictApp[node[0]]) >= int(numberTips*conf):
            g.node(node[0],color=colors[0])
        elif node[1][0] != '':
            g.node(node[0],color=colors[1])
        else:
            g.node(node[0],color=colors[2])

    for node in s:
            for neib in node[1]:
                if neib != '':
                    g.edge(node[0],neib,dir="back")
    g.render(path +'\im\Tangle' + t + 'sec')
    os.remove(path +'\im\Tangle' + t + 'sec')

os.remove(outfile)