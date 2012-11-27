import sys, string, re

original_filename = sys.argv[1]
no_ext_file_name = original_filename.replace(".txt","")
cleaned_up_filename = no_ext_file_name + "_clean.txt"

f = open(cleaned_up_filename, 'wb')

g = open("rejects.txt", 'wb')

for line in open(original_filename, 'r'):
	line_new = line.replace("[]","0.0")
	s = [m.start() for m in re.finditer(' ', line_new)]
	
	if ( ((len(line_new) == 24) and (s[0] == 2) and (s[1] == 12)) or ((len(line_new) == 25) and (s[0] == 3) and (s[1] == 13)) ):
		f.write(line_new)


