class Write(webapp2.RequestHandler):
    def get(self):

        run_name     = cgi.escape(self.request.get('run_name'))
        date         = cgi.escape(self.request.get('date'))
        time         = cgi.escape(self.request.get('time'))
        location     = cgi.escape(self.request.get('location'))
        heading      = cgi.escape(self.request.get('heading'))
        altitude     = cgi.escape(self.request.get('altitude'))
        rfid_tag     = cgi.escape(self.request.get('rfid_tag'))
        tag_strength = cgi.escape(self.request.get('tag_strength'))

        if user:
            self.response.headers['Content-Type'] = 'text/plain'
            self.response.out.write('Hello, ' + user.nickname())
        else:
            self.redirect(users.create_login_url(self.request.uri))

app = webapp2.WSGIApplication([('/', MainPage)],
                              debug=True)

class Read(webapp2.RequestHandler):
    def get(self):
        user = users.get_current_user()

        if user:
            self.response.headers['Content-Type'] = 'text/plain'
            self.response.out.write('Hello, ' + user.nickname())
        else:
            self.redirect(users.create_login_url(self.request.uri))

app = webapp2.WSGIApplication([('/', MainPage)],
                              debug=True)