#We acknowledge use of the CHIME/FRB Public Database, provided at https://www.chime-frb.ca/ by the CHIME/FRB Collaboration.

#search through the CHIME repeating FRBs list and plot them on a timeline with pygame library
import json
import pygame
from datetime import datetime, timedelta, timezone
import time

class FRB:
    def __init__(self,timestamp,name, id):
        date_object = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
        numbertime = int(date_object.timestamp())
        self.timestamp = numbertime
        self.name = name
        self.id = id
    

count = 1
frbs = []
totalBursts = {}
last2Months = {}
with open("repeaters.json", "r") as read_file:
    data = json.load(read_file)
    for frbName, frb in data.items():
        count+=1
        for eventName,event in frb.items():
            try:
                timestamp = event["timestamp"]["value"]
                thisBurst = FRB(timestamp,frbName,count)
                frbs.append(thisBurst)
                totalBursts[frbName] = totalBursts.get(frbName, 0) + 1
                date_object = datetime.strptime(timestamp, "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
                currenttime = datetime.now(timezone.utc)
                # Calculate the difference
                two_months = timedelta(days=60)  # Approximating 2 months as 60 days
                if (currenttime - date_object < two_months):
                    last2Months[frbName] = last2Months.get(frbName, 0) + 1
            except:
                print("frb event does not have timestamp: "+frbName+" "+eventName)



pygame.init()
pygame.font.init()
pygame.display.set_caption("Repeating FRB Timeline")
font = pygame.font.SysFont('arial',10)
screen = pygame.display.set_mode([1000,500])
screen.fill((255,255,255))
earliest = frbs[0].timestamp
latest = frbs[0].timestamp
for frb in frbs:
    if frb.timestamp < earliest:
        earliest = frb.timestamp
    if frb.timestamp > latest:
        latest = frb.timestamp


scores = {}
for frbName in totalBursts.keys():
    scores[frbName] = totalBursts.get(frbName)
sorted = dict(sorted(scores.items(), key=lambda item: item[1], reverse=True))
print(sorted)

for frb in frbs:
    pygame.draw.circle(screen,(0,0,0),((frb.timestamp-earliest)/(latest-earliest)*900+95,(frb.id/count)*475),1)
    text_surface = font.render(frb.name +"  = "+ str(scores[frb.name]),False,(0,0,0))
    screen.blit(text_surface,(3,(frb.id/count)*475-6))



pygame.display.update()

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

pygame.quit()
