import re

def listToSite(listStr):
	outputString = 'site';
	coords = listStr.split(',');
	digitRegex = '-?[\s]?[0-9]+';
	digits = [0]*len(coords)
	for ii in range(len(coords)):
		digit = re.findall(digitRegex, coords[ii])
		if len(digit) != 0:
			num = int(digit[0]);
			if num > 0:
				sgn = '+'
			else:
				sgn = '-'
			outputString += (sgn + str(ii))*abs(num)
	return outputString

def moveSubscript(subscriptString, bracketString):
	subscriptList = subscriptString.split(',')
	return subscriptList[0] + '(' + bracketString + ', ' + subscriptList[1] + ')'


def mathematicaToC(string):
	substitutions = {
	"Cos" : "cos",
	"Sin" : "sin",
	"Tan" : "tan",
	"Sqrt" : "sqrt",
	"Power" : "pow",
	"Exp" : "exp"
	}
	listRegex = 'List\((.*?)\)';
	string = re.sub(listRegex, lambda x : listToSite(x.group(0)), string)
	subscriptRegex = 'Subscript\((.*?)\)\((.*?)\)'
	string = re.sub(subscriptRegex, lambda x : moveSubscript(x.group(1), x.group(2)), string)
	for mathExpr, cExpr in substitutions.items():
		string = string.replace(mathExpr, cExpr)
	return string + ';'

print(mathematicaToC('''Case 0: 
-Subscript(field,0)(List(x,-1 + y,z)) - Subscript(field,0)(List(x,y,-1 + z)) + 4*Subscript(field,0)(List(x,y,z)) - Subscript(field,0)(List(x,y,1 + z)) - Subscript(field,0)(List(x,1 + y,z)) + Subscript(field,1)(List(x,-1 + y,z)) - Subscript(field,1)(List(x,y,z)) - Subscript(field,1)(List(1 + x,-1 + y,z)) + Subscript(field,1)(List(1 + x,y,z)) + Subscript(field,2)(List(x,y,-1 + z)) - Subscript(field,2)(List(x,y,z)) - Subscript(field,2)(List(1 + x,y,-1 + z)) + Subscript(field,2)(List(1 + x,y,z)) + 2*e*Sin(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,3)(List(x,y,z))*Subscript(field,3)(List(1 + x,y,z)) - 2*e*Cos(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,3)(List(1 + x,y,z))*Subscript(field,4)(List(x,y,z)) + 2*e*Cos(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,3)(List(x,y,z))*Subscript(field,4)(List(1 + x,y,z)) + 2*e*Sin(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,4)(List(x,y,z))*Subscript(field,4)(List(1 + x,y,z))
Case 1: 
Subscript(field,0)(List(-1 + x,y,z)) - Subscript(field,0)(List(-1 + x,1 + y,z)) - Subscript(field,0)(List(x,y,z)) + Subscript(field,0)(List(x,1 + y,z)) - Subscript(field,1)(List(-1 + x,y,z)) - Subscript(field,1)(List(x,y,-1 + z)) + 4*Subscript(field,1)(List(x,y,z)) - Subscript(field,1)(List(x,y,1 + z)) - Subscript(field,1)(List(1 + x,y,z)) + Subscript(field,2)(List(x,y,-1 + z)) - Subscript(field,2)(List(x,y,z)) - Subscript(field,2)(List(x,1 + y,-1 + z)) + Subscript(field,2)(List(x,1 + y,z)) + 2*e*Sin(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,3)(List(x,y,z))*Subscript(field,3)(List(x,1 + y,z)) - 2*e*Cos(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,3)(List(x,1 + y,z))*Subscript(field,4)(List(x,y,z)) + 2*e*Cos(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,3)(List(x,y,z))*Subscript(field,4)(List(x,1 + y,z)) + 2*e*Sin(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,4)(List(x,y,z))*Subscript(field,4)(List(x,1 + y,z))
Case 2: 
Subscript(field,0)(List(-1 + x,y,z)) - Subscript(field,0)(List(-1 + x,y,1 + z)) - Subscript(field,0)(List(x,y,z)) + Subscript(field,0)(List(x,y,1 + z)) + Subscript(field,1)(List(x,-1 + y,z)) - Subscript(field,1)(List(x,-1 + y,1 + z)) - Subscript(field,1)(List(x,y,z)) + Subscript(field,1)(List(x,y,1 + z)) - Subscript(field,2)(List(-1 + x,y,z)) - Subscript(field,2)(List(x,-1 + y,z)) + 4*Subscript(field,2)(List(x,y,z)) - Subscript(field,2)(List(x,1 + y,z)) - Subscript(field,2)(List(1 + x,y,z)) + 2*e*Sin(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,3)(List(x,y,z))*Subscript(field,3)(List(x,y,1 + z)) - 2*e*Cos(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,3)(List(x,y,1 + z))*Subscript(field,4)(List(x,y,z)) + 2*e*Cos(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,3)(List(x,y,z))*Subscript(field,4)(List(x,y,1 + z)) + 2*e*Sin(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,4)(List(x,y,z))*Subscript(field,4)(List(x,y,1 + z))
Case 3: 
-2*(Cos(e*Subscript(field,0)(List(-1 + x,y,z)))*Subscript(field,3)(List(-1 + x,y,z)) + Cos(e*Subscript(field,1)(List(x,-1 + y,z)))*Subscript(field,3)(List(x,-1 + y,z)) + Cos(e*Subscript(field,2)(List(x,y,-1 + z)))*Subscript(field,3)(List(x,y,-1 + z)) - 6*Subscript(field,3)(List(x,y,z)) + 2*selfCoupling*Power(vev,2)*Subscript(field,3)(List(x,y,z)) - 2*selfCoupling*Power(Subscript(field,3)(List(x,y,z)),3) + Cos(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,3)(List(x,y,1 + z)) + Cos(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,3)(List(x,1 + y,z)) + Cos(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,3)(List(1 + x,y,z)) + Sin(e*Subscript(field,0)(List(-1 + x,y,z)))*Subscript(field,4)(List(-1 + x,y,z)) + Sin(e*Subscript(field,1)(List(x,-1 + y,z)))*Subscript(field,4)(List(x,-1 + y,z)) + Sin(e*Subscript(field,2)(List(x,y,-1 + z)))*Subscript(field,4)(List(x,y,-1 + z)) - 2*selfCoupling*Subscript(field,3)(List(x,y,z))*Power(Subscript(field,4)(List(x,y,z)),2) - Sin(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,4)(List(x,y,1 + z)) - Sin(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,4)(List(x,1 + y,z)) - Sin(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,4)(List(1 + x,y,z)))
Case 4: 
2*(Sin(e*Subscript(field,0)(List(-1 + x,y,z)))*Subscript(field,3)(List(-1 + x,y,z)) + Sin(e*Subscript(field,1)(List(x,-1 + y,z)))*Subscript(field,3)(List(x,-1 + y,z)) + Sin(e*Subscript(field,2)(List(x,y,-1 + z)))*Subscript(field,3)(List(x,y,-1 + z)) - Sin(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,3)(List(x,y,1 + z)) - Sin(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,3)(List(x,1 + y,z)) - Sin(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,3)(List(1 + x,y,z)) - Cos(e*Subscript(field,0)(List(-1 + x,y,z)))*Subscript(field,4)(List(-1 + x,y,z)) - Cos(e*Subscript(field,1)(List(x,-1 + y,z)))*Subscript(field,4)(List(x,-1 + y,z)) - Cos(e*Subscript(field,2)(List(x,y,-1 + z)))*Subscript(field,4)(List(x,y,-1 + z)) + 6*Subscript(field,4)(List(x,y,z)) - 2*selfCoupling*Power(vev,2)*Subscript(field,4)(List(x,y,z)) + 2*selfCoupling*Power(Subscript(field,3)(List(x,y,z)),2)*Subscript(field,4)(List(x,y,z)) + 2*selfCoupling*Power(Subscript(field,4)(List(x,y,z)),3) - Cos(e*Subscript(field,2)(List(x,y,z)))*Subscript(field,4)(List(x,y,1 + z)) - Cos(e*Subscript(field,1)(List(x,y,z)))*Subscript(field,4)(List(x,1 + y,z)) - Cos(e*Subscript(field,0)(List(x,y,z)))*Subscript(field,4)(List(1 + x,y,z)))'''))