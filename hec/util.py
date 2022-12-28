import sys

def printexit(exitmessage):
    print(exitmessage)
    sys.exit()


def wraptotext(textinput, size=None):
    if size is None:
        size = 120
    textlist = wrap(textinput, size)
    textresult = textlist[0]
    for itext in range(1, len(textlist)):
        textresult += '\n' + textlist[itext]
    return textresult


def timenow():
    now = datetime.now()
    return now.strftime("%m/%d/%Y, %H:%M:%S") + " UTC"


def float32fromstrwithNaN(instr):
    if instr == 'NaN':
        return NaN
    return np.float32(instr)


def strrnd(value):
    return str(round(value, 4))
