# check that class registration works
try: import oofem
except: import oofempy as oofem

def test_0():
    kk=oofem.getClassFactory().getRegisteredNames()
    em=kk['EngngModel']
    assert 'linearstatic' in em
    assert 'nonlinearstatic' in em
    assert 'nonstationaryproblem' in em

    try: from rich.pretty import pprint
    except: from pprint import pprint
    pprint(kk)
if __name__=='__main__':
    test_0()

