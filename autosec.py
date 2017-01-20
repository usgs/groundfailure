#!/usr/bin/env python

#stdlib imports
import configparser
import os.path
import urllib.request, urllib.error, urllib.parse
import json
import smtplib
import sys
import sqlite3
from collections import OrderedDict
from datetime import datetime
import email
from email import encoders
from email.mime.text import MIMEText
#from email.message import Message
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
#from email.mime.application import MIMEApplication
from email.mime.base import MIMEBase
import mimetypes

#third party
from impactutils.io.cmd import get_command_output

CONFIGFILE = 'mailconfig.ini'
CONFIGLIST = 'configlist.txt'
#dictionary containing table definitions and column definitions column:type
tablecols = [('id', 'integer primary key'),
             ('eventcode', 'text'),
             ('version', 'integer'),
             ('lat', 'real'),
             ('lon', 'real'),
             ('depth', 'real'),
             ('time', 'timestamp'),
             ('mag', 'real'),
             ('alert', 'text'),
             ('maxmmi', 'real'),
             ('location', 'text')]
TABLES = {'shakemap': OrderedDict(tablecols)}
DBFILE = 'mail.db'

FEED = 'http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_hour.geojson'
#FEED = 'http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_day.geojson'
#FEED = 'http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_week.geojson'

ALERTLEVELS = ['green', 'yellow', 'orange', 'red', 'pending']


def mailUsers(filenames, event, config, filetypes=('.png', '.html')):
    eid = event['eventcode']
    title = event['title']
    vnum = event['version']
    filepath = os.path.dirname(filenames[0])
    text = """Attached are the most recent secondary hazard files for ShakeMap v%i of event id %s.
    Event name: %s
    If html file is attached, view it in any web browser.

    Only some of the output files are sent through email, the rest of the files, including GIS files and
    pdfs can be found here: %s

    Do not reply to this message.

    Contact kallstadt@usgs.gov with questions.""" % (vnum, eid, title, filepath)
    subject = 'Groundfailure Maps for v%i of %s' % (vnum, title)
    server = config.get('MAIL', 'server')
    sender = config.get('MAIL', 'sender')
    recipients = config.get('MAIL', 'recipients').split(',')
    session = smtplib.SMTP(server)
    for recipient in recipients:
        outer = MIMEMultipart()
        outer['Subject'] = subject
        outer['To'] = recipient
        outer['From'] = sender
        outer['Date'] = email.utils.formatdate()
        outer.attach(MIMEText(text))
        for filen in filenames:
            if os.path.splitext(filen)[1] in filetypes:  # only attach specified extensions
                ctype, encoding = mimetypes.guess_type(filen)
                if ctype is None or encoding is not None:
                    # No guess could be made, or the file is encoded (compressed), so
                    # use a generic bag-of-bits type.
                    ctype = 'application/octet-stream'
                maintype, subtype = ctype.split('/', 1)
                if maintype == 'text':
                    fp = open(filen)
                    # Note: we should handle calculating the charset
                    msg = MIMEText(fp.read(), _subtype=subtype)
                    fp.close()
                elif maintype == 'image':
                    fp = open(filen, 'rb')
                    msg = MIMEImage(fp.read(), _subtype=subtype)
                    fp.close()
                else:
                    fp = open(filen, 'rb')
                    msg = MIMEBase(maintype, subtype)
                    msg.set_payload(fp.read())
                    fp.close()
                    # Encode the payload using Base64
                    encoders.encode_base64(msg)
                # Set the filename parameter
                msg.add_header('Content-Disposition', 'attachment', filename=os.path.basename(filen))
                outer.attach(msg)

        msgtxt = outer.as_string()
        session.sendmail(sender, recipient, msgtxt)

    session.quit()


def runGF(modelconfig, shakefile):
    cmd = 'gfail --gis -pn -pi %s %s' % (modelconfig, shakefile)
    retcode, stdout, stderr = get_command_output(cmd)
    temp = stdout.decode('utf-8')
    if temp.find('Files created:\n') > -1:
        temp = temp.split('Files created:\n')[1]
        filenames = []
        for line in temp.split('\n'):
            if 'pdf' in line or 'png' in line or 'html' in line or 'hdf5' in line or 'bil' in line:
                filenames.append(line)
    else:
        raise Exception('Did not find any files output by runGF, these warnings were output: \
                        %s\n' % stderr)
        return
    return filenames, stdout


def getProductInfo(shakemap, pager):
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
    fh = urllib.request.urlopen(FEED)
    data = fh.read().decode('utf8')
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
        fh = urllib.request.urlopen(detailurl)
        data = fh.read().decode('utf8')
        jdict2 = json.loads(data)
        fh.close()
        shakemap = jdict2['properties']['products']['shakemap'][0]
        pager = jdict2['properties']['products']['losspager'][0]
        #check pager thresholds
        pmag = float(pager['properties']['magnitude'])
        pmmi = float(pager['properties']['maxmmi'])
        palert = ALERTLEVELS.index(pager['properties']['alertlevel'])
        getShake = False
        if 'mag' in thresholds and pmag > thresholds['mag']:
            getShake = True
        if 'mmi' in thresholds and pmmi > thresholds['mmi']:
            getShake = True
        if 'eis' in thresholds and palert >= ALERTLEVELS.index(thresholds['eis']):
            getShake = True
        if getShake:
            edict = getProductInfo(shakemap, pager)
            edict['time'] = datetime.utcfromtimestamp(event['properties']['time']/1000)
            edict['title'] = event['properties']['title']
            eventlist.append(edict.copy())
    return eventlist


def connect():
    dbfile = os.path.join(os.path.expanduser('~'), DBFILE)
    doCreate = False
    if not os.path.isfile(dbfile):
        doCreate = True

    db = sqlite3.connect(dbfile)
    cursor = db.cursor()
    if doCreate:
        for tablename, table in TABLES.items():
            querynuggets = []
            for key, value in table.items():
                nugget = '%s %s' % (key, value)
                querynuggets.append(nugget)
            query = 'CREATE TABLE %s (%s)' % (tablename, ', '.join(querynuggets))
            cursor.execute(query)
            db.commit()

    return (db, cursor)


def main():
    print('%s - Running autosec' % datetime.now())
    configfile = os.path.join(os.path.expanduser('~'), CONFIGFILE)
    config = configparser.ConfigParser()
    config.read_file(open(configfile))
    sections = config.sections()
    if 'THRESHOLDS' not in sections or 'MAIL' not in sections:
        print('Missing THRESHOLDS or MAIL section in %s.  Returning' % configfile)
        sys.exit(1)
    thresh = {}
    for key in config.options('THRESHOLDS'):
        value = config.get('THRESHOLDS', key)
        if key in ['mag', 'mmi']:
            value = float(value)
        thresh[key] = value

    db, cursor = connect()
    recentevents = getRecentEvents(thresh)
    for event in recentevents:
        fmt = 'SELECT id FROM shakemap WHERE eventcode="%s" AND version=%i AND time="%s"'
        query = fmt % (event['eventcode'], event['version'], event['time'])
        cursor.execute(query)
        row = cursor.fetchone()
        if row is None:
            #this event has not been processed before
            modelconfig = os.path.join(os.path.expanduser('~'), CONFIGLIST)
            filenames, stdout = runGF(modelconfig, event['url'])
            print(stdout.decode('utf-8'))
            if not filenames:
                print('No outputs found, problem with codes\n')
                continue
            mailUsers(filenames, event, config, filetypes=('.png', '.html'))
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
            insertquery = fmt % (eid, enum, elat, elon, edepth, str(etime), emag, alert, maxmmi, eloc)
            cursor.execute(insertquery)
            db.commit()

if __name__ == '__main__':
    main()
