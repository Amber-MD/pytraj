import uuid
from IPython.display import HTML, Javascript, display

pb = """
<div id="{id0}">{msg}</div>
<div style="border: 1px solid black; width:500px">
  <div id="{id1}" style="background-color:{color}; width:0%">&nbsp;</div>
</div>
"""


def init_display(divid0, divid1, color='#0080FF'):
    my_html = pb.format(id0=divid0, id1=divid1, color=color, msg='0')
    display(HTML(my_html))

def make_bar(idx, max_frames, divid0, divid1):
    '''work with jupyter notebook.
    '''
    percent = str(100*idx/max_frames)
    s0 = "$('div#{id0}').text({percent});".format(id0=divid0, percent=percent)
    s1 = "$('div#{id1}').width('{percent}%');".format(id1=divid1, percent=percent)
    s = s0 + '\n' + s1
    display(Javascript(s))

def log_progress(sequence, every=1, size=None, color='#0080FF'):
    # adapt from: https://github.com/alexanderkuk/log-progress
    from IPython.display import display


    index = 0
    divid0 = str(uuid.uuid4())
    divid1 = str(uuid.uuid4())
    init_display(divid0, divid1, color=color)

    for index, record in enumerate(sequence, 1):
        if index == 1 or index % every == 0:
            make_bar(index, size, divid0, divid1)
        yield record


class ProgressBarTrajectory(object):
    def __init__(self, traj, every=1, color="#0080FF"):
        self.traj = traj
        self.every = every
        self.color = color

        for att in dir(traj):
            if not att.startswith('__'):
                setattr(self, att, getattr(traj, att))
            if att in ['__getstate__', '__setstate__', '_split_iterators']:
                setattr(self, att, getattr(traj, att))

    def __getitem__(self, index):
        return self.traj[index]

    def __iter__(self):
        for frame in log_progress(self.traj, every=self.every,
                size=self.n_frames, color=self.color):
            yield frame
