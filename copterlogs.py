import webapp2, cgi
from google.appengine.api import users
from google.appengine.ext import db

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
        # user = users.get_current_user()

        # if user:
        #     copterlog.put()
        # else:
        #     self.redirect(users.create_login_url(self.request.uri))

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