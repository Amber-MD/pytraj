import uuid
from IPython.display import HTML, Javascript, display

divid = str(uuid.uuid4())

pb = HTML(
"""
<div style="border: 1px solid black; width:500px">
  <div id="%s" style="background-color:blue; width:0%%">&nbsp;</div>
</div> 
""" % divid)

def init_display():
    display(pb)

def make_bar(idx, max_frames):
    '''work with jupyter notebook.
    '''
    s = "$('div#%s').width('%i%%')" % (divid, 100*idx/max_frames)
    display(Javascript(s))
