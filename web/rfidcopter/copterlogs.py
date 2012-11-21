import webapp2, cgi
from google.appengine.api import users
from google.appengine.ext import db
import simplejson

class CopterLog(db.Model):
    run_name     = db.StringProperty() # name of run, if any
    date         = db.DateProperty(auto_now_add=True) # Date of log
    time         = db.TimeProperty(auto_now_add=True) # Time of log
    location     = db.GeoPtProperty() # Lat/Long of copter
    heading      = db.IntegerProperty() # 0-359, degrees from North
    altitude     = db.FloatProperty() # altitude, in ft, imperial system ftw
    rfid_tag     = db.StringProperty() # RFID UID
    tag_strength = db.FloatProperty() # Strength of the tag

def copterlog_key(copterlog_name=None):
    """Constructs a Datastore key for a Guestbook entity with guestbook_name."""
    return db.Key.from_path('CopterLog', copterlog_name or 'default_copterlog')

class Write(webapp2.RequestHandler):
    def get(self):
        copterlog = CopterLog(parent=copterlog_key())
        copterlog.run_name     = cgi.escape(self.request.get('run_name'))
        # copterlog.date         = cgi.escape(self.request.get('date'))
        # copterlog.time         = cgi.escape(self.request.get('time'))
        copterlog.location     = cgi.escape(self.request.get('location'))
        copterlog.heading      = int(cgi.escape(self.request.get('heading')))
        copterlog.altitude     = float(cgi.escape(self.request.get('altitude')))
        copterlog.rfid_tag     = cgi.escape(self.request.get('rfid_tag'))
        copterlog.tag_strength = float(cgi.escape(self.request.get('tag_strength')))
        copterlog.put() 

    def post(self):
        json_obj = simplejson.loads(self.request.body)

        for point in json_obj:
            copterlog = CopterLog(parent=copterlog_key())
            copterlog.run_name     = point["run_name"]
            # copterlog.date         = point["date"]
            # copterlog.time         = point["time"]
            copterlog.location     = point["location"]
            copterlog.heading      = int(point["heading"])
            copterlog.altitude     = float(point["altitude"])
            copterlog.rfid_tag     = point["rfid_tag"]
            copterlog.tag_strength = float(point["tag_strength"])
            copterlog.put()

        # print f
        # print f

def gql_json_parser(query_obj):
    result = []
    for entry in query_obj:
        result.append(dict([(p, unicode(getattr(entry, p))) for p in entry.properties()]))
    return result

class Read(webapp2.RequestHandler):
    def get(self):
        run_name = cgi.escape(self.request.get('run_name'))
        logs = CopterLog.all()

        if run_name:
            logs.filter("run_name =", run_name)
        
        json_query_data =  gql_json_parser(logs)
        self.response.headers['Content-Type'] = 'application/json'
        self.response.out.write(simplejson.dumps(json_query_data))


app = webapp2.WSGIApplication([
                ('/write/', Write),
                ('/read/', Read)
                ], debug=True)