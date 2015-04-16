import os
import sys


def autop(x):
    print x
    a = os.popen('autopep8 ' + x + ' -d').read()
    if a == '':
        print '>>>> No change.'
    else:
        print a
        print x
        a = raw_input("Do you want to change?(y/n):")
        if a == 'y':
            a = os.popen('autopep8 ' + x + ' -i').read()
        else:
            print "Didn't change it."


def check(x):
    l = os.listdir(x)
    for i in l:
        if os.path.isfile(os.path.join(x, i)) and os.path.join(x, i).endswith('.py') and not i.startswith('.'):
            autop(os.path.join(x, i))
        if os.path.isdir(os.path.join(x, i)) and not i.startswith('.') and i != 'build':
            check(os.path.join(x, i))

print '#' * int(os.popen('stty size').read().split()[-1])
check(os.path.abspath(os.path.dirname(sys.argv[0])))
print '#' * int(os.popen('stty size').read().split()[-1])
