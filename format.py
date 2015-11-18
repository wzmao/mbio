import os
import sys


def autop(x, allyes=0):
    print x
    a = os.popen('autopep8 ' + x + ' -d').read()
    if a == '':
        print '>>>> No change.'
    else:
        print a
        print x
        if not allyes:
            a = raw_input("Do you want to change?(y/n):")
            if a == 'y':
                a = os.popen('autopep8 ' + x + ' -i').read()
            else:
                print "Didn't change it."
        else:
            print 'Allyes=1 so correct it automatically.'
            a = os.popen('autopep8 ' + x + ' -i').read()


def check(x, allyes=0):
    l = os.listdir(x)
    for i in l:
        if os.path.isfile(os.path.join(x, i)) and os.path.join(x, i).endswith('.py') and not i.startswith('.'):
            autop(os.path.join(x, i), allyes=allyes)
        if os.path.isdir(os.path.join(x, i)) and not i.startswith('.') and i != 'build':
            check(os.path.join(x, i), allyes=allyes)

print '#' * int(os.popen('stty size').read().split()[-1])
if len(sys.argv) > 1 and any([i in sys.argv[1:] for i in ['y', 'Y', '-y', '-Y']]):
    allyes = 1
else:
    allyes = 0
check(os.path.abspath(os.path.dirname(sys.argv[0])), allyes=allyes)
print '#' * int(os.popen('stty size').read().split()[-1])
