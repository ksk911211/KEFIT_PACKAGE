def read_namelist_str(line,var_name,vars,vartype):

	line = line.split('!')[0]
	line = line.split('#')[0]
	line = line.split('\n')[0]

	try:
		name = line.split('=')[0].split()[0].lower()
	except:
		name = None
	
	if (name == var_name.lower()):

		if (vartype < 4):
			if (vartype == 1):
				var = int(float(line.split('=')[1].replace(' ','')))
			elif (vartype == 2):
				var = float(line.split('=')[1].replace(' ',''))
			elif (vartype == 3):
				var = line.split('=')[1].replace(' ','')
				
			exist = 1
		else:
			var2 = line.split('=')[1].replace(' ','').lower()
			
			if (var2 == 'true'):
				var = True
				
			elif (var2 == 'false'):
				var = False
			else:
				var = vars
	else:
		var = vars;

	line2 = line.split('!')[0].split()
	if(line2 == []):
		var = vars;

	if (var == ''):
		var = None
	
	return (var)
