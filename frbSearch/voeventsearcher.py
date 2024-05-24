#We acknowledge use of the CHIME/FRB Public Database, provided at https://www.chime-frb.ca/ by the CHIME/FRB Collaboration.

import json
from dateutil import parser #pip install python-dateutil
import time

#search through the chime VOEvent database .json file for FRBs
#matching a criteria defined in the variables beneath this
#and then print out the results
dec = 34
decRange = 2

ra = 180
raRange = 180

dm = 500
dmRange = 500

timeWindow = [int(parser.parse("2024-02-01").timestamp()),int(parser.parse("2024-05-24").timestamp())]

class FRB:
    def __init__(self,event_id,alert_type,detected_time,dm,snr,ra,dec,localization_error):
        self.event_id = event_id
        self.alert_type = alert_type
        
        date_object = parser.parse(detected_time)
        numbertime = int(date_object.timestamp())
        self.detected_time = numbertime
        
        self.dm = float(dm)
        self.snr = float(snr)
        self.ra = float(ra)
        self.dec = float(dec)
        self.localization_error = float(localization_error)
    
frbs = []

#load data
with open("chimefrb_voevent_data.json", "r") as read_file:
    data = json.load(read_file)
    for record in data:
        for event in record["records"]:
            try:
                burst = FRB(
                    record["event_id"],
                    event["Alert_Type"],
                    event["Detected"],
                    event["DM"],
                    event["SNR"],
                    event["RA"],
                    event["Dec"],
                    event["Localization_Error"]
                    )
                frbs.append(burst)
            except:
                print("data not complete")
            
def validate(f):
    
    if f.dec < dec-decRange or f.dec > dec+decRange:
        print("out of DEC range: " + str(f.dec))
        return False

    if f.ra < ra-raRange or f.ra > ra+raRange:
        print("out of RA range: " + str(f.ra))
        return False
    
    if f.dm < dm-dmRange or f.dm > dm+dmRange:
        print("out of DM range: " + str(f.dm))
        return False
    
    
    if f.detected_time < timeWindow[0] or f.detected_time > timeWindow[1]:
        print("out of time window: " + str(f.detected_time))
        return False
    
    return True


candidates = []

for f in frbs:
    if validate(f):
        candidates.append(f)

for c in candidates:
    print("found candidate: " + c.event_id)