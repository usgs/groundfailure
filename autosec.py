#!/usr/bin/env python

#stdlib imports
import ConfigParser
import os.path
import urllib2
import json
import smtplib
import sys
import sqlite3
from collections import OrderedDict
from datetime import datetime
import email
from email import encoders
from email.mime.text import MIMEText
from email.message import Message
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.mime.base import MIMEBase

#third party
from neicio.cmdoutput import getCommandOutput

CONFIGFILE = 'mailconfig.ini'
#dictionary containing table definitions and column definitions column:type
tablecols  = [('id','integer primary key'),
              ('eventcode','text'),
              ('version','integer'),
              ('lat','real'),
              ('lon','real'),
              ('depth','real'),
              ('time','timestamp'),
              ('mag','real'),
              ('alert','text'),
              ('maxmmi','real'),
              ('location','text')]
TABLES = {'shakemap':OrderedDict(tablecols)}
DBFILE = 'mail.db'

FEED = 'http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_hour.geojson'
#FEED = 'http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_day.geojson'
#FEED = 'http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_week.geojson'

ALERTLEVELS = ['green','yellow','orange','red','pending']

SECHAZ = 'sechaz.py'

def mailUsers(pdf,png,event,config):
    eid = event['eventcode']
    vnum = event['version']
    text = 'Attached are the two most recent secondary hazard pdf/png files for v%i of %s.' % (vnum,eid)
    subject = 'Secondary Hazard Maps for v%i of %s' % (vnum,eid)
    server = config.get('MAIL','server')
    sender = config.get('MAIL','sender')
    recipients  = config.get('MAIL','recipients').split(',')
    session = smtplib.SMTP(server)
    for recipient in recipients:
        outer = MIMEMultipart()
        outer['Subject'] = subject
        outer['To'] = recipient
        outer['From'] = sender
        outer['Date'] = email.utils.formatdate()
        firstSubMsg=Message()
        firstSubMsg["Content-type"]="text/plain"
        firstSubMsg["Content-transfer-encoding"]="7bit"
        firstSubMsg.set_payload(text)
        outer.attach(firstSubMsg)
        fp = open(png, 'rb')
        pngmsg = MIMEImage(fp.read(), _subtype=None)
        fp.close()
        pngmsg.add_header('Content-Disposition', 'attachment', filename=os.path.basename(png))
        outer.attach(pngmsg)

        fp = open(pdf, 'rb')
        pdfmsg = MIMEBase('application/pdf', None)
        pdfmsg.set_payload(fp.read())
        fp.close()
        pdfmsg.add_header('Content-Disposition', 'attachment', filename=os.path.basename(pdf))
        #Encode the payload using Base64
        encoders.encode_base64(pdfmsg)
        
        outer.attach(pdfmsg)
        msgtxt = outer.as_string()
        session.sendmail(sender,recipient, msgtxt)
    
    session.quit()
        
def runSecondary(url,thisdir):
    cmd = os.path.join(thisdir,SECHAZ)
    cmd += ' %s' % url
    retcode,stdout,stderr = getCommandOutput(cmd)
    pdf = None
    png = None
    for line in stdout.split('\n'):
        if line.find('Saving map output') > -1:
            parts = line.split()
            tfile = parts[-1]
            if tfile.find('pdf') > -1:
                pdf = tfile
            if tfile.find('png') > -1:
                png = tfile
    return (pdf,png)

def getProductInfo(shakemap,pager):
    edict = {}
    edict['eventcode'] = shakemap['code']
    edict['version'] = int(shakemap['properties']['version'])
    edict['lat'] = float(shakemap['properties']['latitude'])
    edict['lon'] = float(shakemap['properties']['longitude'])
    edict['depth'] = float(shakemap['properties']['depth'])
    #shakemap properties don't have event time, fill this in from event info
    edict['mag'] = float(shakemap['properties']['magnitude'])
    edict['location'] = shakemap['properties']['event-description']
    edict['url'] = shakemap['contents']['download/grid.xml']['url']
    edict['alert'] = pager['properties']['alertlevel']
    edict['maxmmi'] = float(pager['properties']['maxmmi'])
    return edict

def getRecentEvents(thresholds):
    fh = urllib2.urlopen(FEED)
    data = fh.read()
    jdict = json.loads(data)
    fh.close()
    eventlist = []
    for event in jdict['features']:
        etypes = event['properties']['types'].strip(',').split(',')
        if 'shakemap' not in etypes:
            continue
        if 'losspager' not in etypes:
            continue
        edict = {}
        detailurl = event['properties']['detail']
        fh = urllib2.urlopen(detailurl)
        data = fh.read()
        jdict2 = json.loads(data)
        fh.close()
        shakemap = jdict2['properties']['products']['shakemap'][0]
        pager = jdict2['properties']['products']['losspager'][0]
        #check pager thresholds
        pmag = float(pager['properties']['magnitude'])
        pmmi = float(pager['properties']['maxmmi'])
        palert = ALERTLEVELS.index(pager['properties']['alertlevel'])
        getShake = False
        if thresholds.has_key('mag') and pmag > thresholds['mag']:
            getShake = True
        if thresholds.has_key('mmi') and pmag > thresholds['mag']:
            getShake = True
        if thresholds.has_key('eis') and palert >= ALERTLEVELS.index(thresholds['eis']):
            getShake = True
        if getShake:
            edict = getProductInfo(shakemap,pager)
            edict['time'] = datetime.utcfromtimestamp(event['properties']['time']/1000)
            eventlist.append(edict.copy())
    return eventlist
        

def connect():
    dbfile = os.path.join(os.path.expanduser('~'),'.secondary',DBFILE)
    doCreate = False
    if not os.path.isfile(dbfile):
        doCreate = True
        
    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    if doCreate:
        for tablename,table in TABLES.iteritems():
            querynuggets = []
            for key,value in table.iteritems():
                nugget = '%s %s' % (key,value)
                querynuggets.append(nugget)
            query = 'CREATE TABLE %s (%s)' % (tablename,', '.join(querynuggets))
            cursor.execute(query)
            db.commit()

    return (db,cursor)

def main():
    print '%s - Running autosec' % datetime.now()
    thisdir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
    configfile = os.path.join(os.path.expanduser('~'),'.secondary',CONFIGFILE)
    config = ConfigParser.ConfigParser()
    config.readfp(open(configfile))
    sections = config.sections()
    if 'THRESHOLDS' not in sections or 'MAIL' not in sections:
        print 'Missing THRESHOLDS or MAIL section in %s.  Returning' % configfile
        sys.exit(1)
    thresh = {}
    for key in config.options('THRESHOLDS'):
        value = config.get('THRESHOLDS',key)
        if key in ['mag','mmi']:
            value = float(value)
        thresh[key] = value

    db,cursor = connect()
    recentevents = getRecentEvents(thresh)
    for event in recentevents:
        fmt = 'SELECT id FROM shakemap WHERE eventcode="%s" AND version=%i AND time="%s"'
        query = fmt % (event['eventcode'],event['version'],event['time'])
        cursor.execute(query)
        row = cursor.fetchone()
        if row is None:
            #this event has not been processed before
            pdf,png = runSecondary(event['url'],thisdir)
            mailUsers(pdf,png,event,config)
            fmt = 'INSERT INTO shakemap (eventcode,version,lat,lon,depth,time,mag,alert,maxmmi,location) VALUES ("%s",%i,%.4f,%.4f,%.1f,"%s",%.1f,"%s",%.1f,"%s")'
            eid = event['eventcode']
            enum = event['version']
            elat = event['lat']
            elon = event['lon']
            edepth = event['depth']
            etime = event['time']
            emag = event['mag']
            alert = event['alert']
            maxmmi = event['maxmmi']
            eloc = event['location']
            insertquery = fmt % (eid,enum,elat,elon,edepth,str(etime),emag,alert,maxmmi,eloc)
            cursor.execute(insertquery)
            db.commit()

if __name__ == '__main__':
    main()
        
