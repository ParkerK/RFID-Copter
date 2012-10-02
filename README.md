RFID-Copter
===========

Reynolds RFID Quadcopter project

Local Install
=============
* Download and install the Google App Engine SDK for Python: https://developers.google.com/appengine/downloads#Google_App_Engine_SDK_for_Python
* Pull this repo, and set up as a new GAE project

Useage
======
* Run the app locally in the GAE Launcher program
* http://localhost:8080/read/ will return all points as a JSON object. Optionally, you can append with a `?run_name=X` to filter the results
* http://localhost:8080/write/ will take a JSON object (sent via POST) and add it to the database. 

Datapoint Properties
====================
run_name     = db.StringProperty() # name of run, if any
date         = db.DateProperty(auto_now_add=True) # Date of log
time         = db.TimeProperty(auto_now_add=True) # Time of log
location     = db.GeoPtProperty() # Lat/Long of copter
heading      = db.IntegerProperty() # 0-359, degrees from North
altitude     = db.FloatProperty() # altitude, in ft, imperial system ftw
rfid_tag     = db.StringProperty() # RFID UID
tag_strength = db.FloatProperty() # Strength of the tag