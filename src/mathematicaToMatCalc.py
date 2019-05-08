# -*- coding: utf-8 -*-
import re

'''
Change prime to apostrophe
Replace minus signs
Remove commas
Remove curly braces
Remove square braces after u or p
Decapitalise Tr
Other square brackets to round
Dots/spaces to asterisks
'''
import sys
print(sys.version)
def mathematicaToMatCalc(string):
	outputString = string;
	outputString = string.replace('^\\[Prime]','\'')
	minusSignRegex = '([,\{])(-)'
	outputString = re.sub(minusSignRegex, lambda x : x.group(1) + 'm', outputString)
	outputString = outputString.replace(',','')
	squareBracketRegex = '([up])\[(.*?)\]'
	outputString = re.sub(squareBracketRegex, lambda x : x.group(1) + x.group(2), outputString)
	outputString = outputString.replace('{','')
	outputString = outputString.replace('}','')
	outputString = outputString.replace('Tr','tr')
	outputString = outputString.replace('[','(')
	outputString = outputString.replace(']',')')
	outputString = outputString.replace('.','*')
	outputString = outputString.replace(' ','*')
	outputString = outputString.replace('u', 'U')
	outputString = outputString.replace('p', 'P')
	return outputString

def matCalcToMathematica(string):
	outputString = string
	outputString = re.sub(u'([UP][0-2m]*)[\'\u22A4]', lambda x : 'ConjugateTranspose[' + x.group(1) + ']', outputString)
	# Awful hack to get compound expressions in transpose; doesn't work for nested brackets
	outputStringReversed = re.sub(u'[\'\u22A4]\)(.*?)\(', lambda x : ']' + x.group(1) + '[esopsnarTetagujnoC', outputString[::-1])
	outputString = outputStringReversed[::-1]
	uGroups = re.findall('U([0-2m]*)', string)
	pCoords = re.findall('P([0-2m]*)', string)
	uCoords = []
	uDirs = []
	for uInfo in uGroups:
		uCoordVec = []
		ii = 0
		while ii < len(uInfo) - 1:
			if uInfo[ii] == 'm':
				uCoordVec.append(-int(uInfo[ii+1]))
				ii = ii+2
			else:
				uCoordVec.append(int(uInfo[ii]))
				ii = ii+1
		uDir = int(uInfo[-1])
		outputString = outputString.replace('U'+uInfo, 'u[{'+','.join(map(str,uCoordVec)) +'},' +str(uDir)+']')
	for pCoord in pCoords:
		pCoordVec = []
		ii = 0
		while ii < len(pCoord):
			if pCoord[ii] == 'm':
				pCoordVec.append(-int(pCoord[ii+1]))
				ii = ii+2
			else:
				pCoordVec.append(int(pCoord[ii]))
				ii = ii+1
		outputString = outputString.replace('P'+pCoord, 'p[{'+','.join(map(str,pCoordVec)) +'}]')
	outputString = re.sub(u'[\*\u22C5\s]', '.', outputString)
	return outputString

def matCalcToLatex(string):
	outputString = string
	uGroups = re.findall('U([0-2m]*)', string)
	pCoords = re.findall('P([0-2m]*)', string)
	uCoords = []
	uDirs = []
	for uInfo in uGroups:
		uCoordVec = []
		ii = 0
		while ii < len(uInfo) - 1:
			if uInfo[ii] == 'm':
				uCoordVec.append(-int(uInfo[ii+1]))
				ii = ii+2
			else:
				uCoordVec.append(int(uInfo[ii]))
				ii = ii+1
		uDir = int(uInfo[-1])
		if uDir == 0:
			dirString = '\\mu'
		else:
			dirString = '\\nu'
		coordString = 'x'
		if uCoordVec[0] == 1:
			coordString += ' + \\mu'
		elif uCoordVec[0] == -1:
			coordString += ' - \\mu'
		if uCoordVec[1] == 1 or uCoordVec[2] == 1:
			coordString += ' + \\nu'
		elif uCoordVec[1] == -1 or uCoordVec[2] == -1:
			coordString += ' + \\nu'
		outputString = outputString.replace('U' + uInfo, 'U_{' + dirString + '}(' + coordString + ')')

	for pCoord in pCoords:
		pCoordVec = []
		ii = 0
		while ii < len(pCoords) - 1:
			if uInfo[ii] == 'm':
				pCoordVec.append(-int(uInfo[ii+1]))
				ii = ii+2
			else:
				pCoordVec.append(int(uInfo[ii]))
				ii = ii+1
		if pCoordVec[0] == 1:
			coordString += ' + \\mu'
		elif pCoordVec[0] == -1:
			coordString += ' - \\mu'
		if pCoordVec[1] == 1 or pCoordVec[2] == 1:
			coordString += ' + \\nu'
		elif pCoordVec[1] == -1 or pCoordVec[2] == -1:
			coordString += ' + \\nu'
		outputString = outputString.replace('P'+pCoord, '\\Phi(' + coordString + ')')
	outputString = re.sub('\u22A4', '^\\dagger ', outputString)
	outputString = re.sub('[\*\u22C5\s]', '', outputString)
	return outputString

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

print(matCalcToMathematica('''
4⋅U0m101⊤⋅U0m100⋅P1m10⊤⋅U0m100⊤⋅(U0m100⋅P1m10⊤⋅U0m100⊤+(−P0m10)⊤+U0m100⋅P1m10⋅U0m100⊤−P0m10)⋅U0m100⋅U1m101+4⋅U0m101⊤⋅(U0m101⋅P000⋅U0m101⊤−P0m10+U0m101⋅P000⊤⋅U0m101⊤+(−P0m10)⊤)⋅U0m101⋅P000⋅U0m101⊤⋅U0m100⋅U1m101+4⋅U00m12⊤⋅U00m10⋅P10m1⊤⋅U00m10⊤⋅(U00m10⋅P10m1⊤⋅U00m10⊤+(−P00m1)⊤+U00m10⋅P10m1⋅U00m10⊤−P00m1)⋅U00m10⋅U10m12+4⋅U00m12⊤⋅(U00m12⋅P000⋅U00m12⊤−P00m1+U00m12⋅P000⊤⋅U00m12⊤+(−P00m1)⊤)⋅U00m12⋅P000⋅U00m12⊤⋅U00m10⋅U10m12+4⋅U0000⋅(U1001⋅U0100⊤⋅U0001⊤+U1002⋅U0010⊤⋅U0002⊤+(U0m100⋅U1m101)⊤⋅U0m101+(U00m10⋅U10m12)⊤⋅U00m12)⋅U0000⋅P100⊤⋅U0000⊤⋅U0000⋅P100⊤+4⋅U0000⋅P100⋅U0000⊤⋅(U0001⋅U0100⋅U1001⊤+U0002⋅U0010⋅U1002⊤+U0m101⊤⋅U0m100⋅U1m101+U00m12⊤⋅U00m10⋅U10m12)⋅U0000⊤⋅U0000⋅P100+4⋅U0000⋅(U1001⋅U0100⊤⋅U0001⊤+U1002⋅U0010⊤⋅U0002⊤+(U0m100⋅U1m101)⊤⋅U0m101+(U00m10⋅U10m12)⊤⋅U00m12)⋅U0000⋅P100⊤⋅U0000⊤⋅U0000⋅P100+4⋅U0000⋅P100⋅U0000⊤⋅(U0001⋅U0100⋅U1001⊤+U0002⋅U0010⋅U1002⊤+U0m101⊤⋅U0m100⋅U1m101+U00m12⊤⋅U00m10⋅U10m12)⋅U0000⊤⋅U0000⋅P100⊤+4⋅(U0000⋅P100⊤⋅U0000⊤+(−P000)⊤+U0000⋅P100⋅U0000⊤−P000)⋅U0000⋅(U1001⋅U0100⊤⋅U0001⊤+U1002⋅U0010⊤⋅U0002⊤+(U0m100⋅U1m101)⊤⋅U0m101+(U00m10⋅U10m12)⊤⋅U00m12)⋅U0000⋅P100⊤+4⋅(U0001⋅U0100⋅U1001⊤+U0002⋅U0010⋅U1002⊤+U0m101⊤⋅U0m100⋅U1m101+U00m12⊤⋅U00m10⋅U10m12)⋅U0000⊤⋅(U0000⋅P100⋅U0000⊤−P000+U0000⋅P100⊤⋅U0000⊤+(−P000)⊤)⋅U0000⋅P100+4⋅(U0000⋅P100⋅U0000⊤−P000+U0000⋅P100⊤⋅U0000⊤+(−P000)⊤)⋅U0000⋅P100⋅U0000⊤⋅(U0001⋅U0100⋅U1001⊤+U0002⋅U0010⋅U1002⊤+U0m101⊤⋅U0m100⋅U1m101+U00m12⊤⋅U00m10⋅U10m12)+4⋅U0001⋅P010⊤⋅U0001⊤⋅(U0001⋅P010⊤⋅U0001⊤+(−P000)⊤+U0001⋅P010⋅U0001⊤−P000)⋅U0001⋅U0100⋅U1001⊤+4⋅U0002⋅P001⊤⋅U0002⊤⋅(U0002⋅P001⊤⋅U0002⊤+(−P000)⊤+U0002⋅P001⋅U0002⊤−P000)⋅U0002⋅U0010⋅U1002⊤+4⋅U0002⋅U0010⋅P101⊤⋅U0010⊤⋅(U0010⋅P101⊤⋅U0010⊤+(−P001)⊤+U0010⋅P101⋅U0010⊤−P001)⋅U0010⋅U1002⊤+4⋅U0001⋅U0100⋅P110⊤⋅U0100⊤⋅(U0100⋅P110⊤⋅U0100⊤+(−P010)⊤+U0100⋅P110⋅U0100⊤−P010)⋅U0100⋅U1001⊤+4⋅U0m101⊤⋅U0m100⋅U1m101⋅P100⊤⋅U1m101⊤⋅(U1m101⋅P100⊤⋅U1m101⊤+(−P1m10)⊤+U1m101⋅P100⋅U1m101⊤−P1m10)⋅U1m101+4⋅U00m12⊤⋅U00m10⋅U10m12⋅P100⊤⋅U10m12⊤⋅(U10m12⋅P100⊤⋅U10m12⊤+(−P10m1)⊤+U10m12⋅P100⋅U10m12⊤−P10m1)⋅U10m12+4⋅U0001⋅U0100⋅U1001⊤⋅(U1001⋅P110⋅U1001⊤−P100+U1001⋅P110⊤⋅U1001⊤+(−P100)⊤)⋅U1001⋅P110⋅U1001⊤+4⋅U0002⋅U0010⋅U1002⊤⋅(U1002⋅P101⋅U1002⊤−P100+U1002⋅P101⊤⋅U1002⊤+(−P100)⊤)⋅U1002⋅P101⋅U1002
	'''))
# matCalcToMathematica('U0001')
