import sys
import textwrap
from datetime import datetime
import numpy as np


startbold = "\033[1m"
resetfonts = "\033[0m"
startred = '\033[31m'

startpurple = '\033[35m'
startyellowbkg = '\033[43m'

def fred(abc):
    print(abc)


NaN = np.float32("NaN")


def printexit(exitmessage):
    print(exitmessage)
    sys.exit()


def wraptotext(textinput, size=None):
    if size is None:
        size = 120
    textlist = textwrap.wrap(textinput, size)
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
