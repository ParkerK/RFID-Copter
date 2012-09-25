from google.appengine.ext import db

class GPS(db.Model):
	run_name = db.StringProperty()
    date = db.DateProperty(required=True) # Date of log
    time = db.TimeProperty(required=True) # Time of log
    location = db.GeoPtProperty(required=True) # Lat/Long of copter
    heading = db.IntegerProperty(required=True) # 0-359°, degrees from North
    altitude = db.FloatProperty(required=True)
    rfid_tag = db.StringProperty(required=True) # RFID UID
    tag_strength = db.FloatProperty(required=True) # Strength of the tag
