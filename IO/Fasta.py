__author__ = 'Wenzhi Mao'
__all__ = ['ReadFasta']

def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()

def ReadFasta(filename, Seq=False):
	f=open(filename)
	fastas=f.read()[1:].split('>')
	f.close()
	name=[fasta.split('\n')[0] for fasta in fastas]
	seq=[''.join(''.join(fasta.split('\n')[1:]).split()) for fasta in fastas]
	if Seq:
		return seq
	else:
		return name , seq


_Startup()