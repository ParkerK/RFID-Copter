import webapp2, cgi
from google.appengine.api import users
from google.appengine.ext import db

class CopterLog(db.Model):
    run_name     = db.StringProperty() # name of run, if any
    date         = db.DateProperty(required=True) # Date of log
    time         = db.TimeProperty(required=True) # Time of log
    location     = db.GeoPtProperty(required=True) # Lat/Long of copter
    heading      = db.IntegerProperty(required=True) # 0-359, degrees from North
    altitude     = db.FloatProperty(required=True) # altitude, in ft, imperial system ftw
    rfid_tag     = db.StringProperty() # RFID UID
    tag_strength = db.FloatProperty() # Strength of the tag

class Write(webapp2.RequestHandler):
    def get(self):
        copterlog = new CopterLog()
        copterlog.run_name     = cgi.escape(self.request.get('run_name'))
        copterlog.date         = cgi.escape(self.request.get('date'))
        copterlog.time         = cgi.escape(self.request.get('time'))
        copterlog.location     = cgi.escape(self.request.get('location'))
        copterlog.heading      = cgi.escape(self.request.get('heading'))
        copterlog.altitude     = cgi.escape(self.request.get('altitude'))
        copterlog.rfid_tag     = cgi.escape(self.request.get('rfid_tag'))
        copterlog.tag_strength = cgi.escape(self.request.get('tag_strength'))

        user = users.get_current_user()

        if user:
            self.response.headers['Content-Type'] = 'text/plain'
            self.response.out.write('Hello, ' + run_name)
            copterlog.put()
        else:
            self.redirect(users.create_login_url(self.request.uri))

class Read(webapp2.RequestHandler):
    def get(self):
        user = users.get_current_user()

        if user:
            self.response.headers['Content-Type'] = 'text/plain'
            self.response.out.write('Hello, ' + user.nickname())
        else:
            self.redirect(users.create_login_url(self.request.uri))

app = webapp2.WSGIApplication([
                ('/write/', Write),
                ('/read/', Read)
                ], debug=True)