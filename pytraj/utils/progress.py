import uuid
from IPython.display import HTML, Javascript, display

divid = str(uuid.uuid4())

pb = """
<div style="border: 1px solid black; width:500px">
  <div id="{0}" style="background-color:{1}; width:0%%">&nbsp;</div>
</div>
"""


def init_display(color='#0080FF'):
    display(HTML(pb.format(divid, color)))


def make_bar(idx, max_frames):
    '''work with jupyter notebook.
    '''
    s = "$('div#%s').width('%i%%')" % (divid, 100 * idx / max_frames)
    display(Javascript(s))
