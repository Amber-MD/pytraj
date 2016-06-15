import os
import uuid
from IPython.display import HTML, Javascript, display

pb = """
<div id="{id0}">{msg}</div>
<div style="border: 1px solid black; width:500px">
  <div id="{id1}" style="background-color:{color}; width:0%">&nbsp;</div>
</div>
"""

circle_html = """
<script>
%s
</script>

<div id="%s"></div>

<script>
    var ele = "#%s";
    $(ele).progressCircle({
        nPercent        : 0,
        showPercentText : true,
        thickness       : 5,
        circleSize      : 60,
    });
</script>
"""

js_circle = """
$("#%s").progressCircle({
    nPercent        : %s,
    showPercentText : true,
    thickness       : 5,
    circleSize      : 60,
});
"""


class CircleProgress(object):
    @classmethod
    def init_display(cls, circle):
        fn = os.path.dirname(__file__) + '/progress-circle/css/circle.css'
        style = "<style>\n" + open(fn).read() + "</style>" 
        fn2 = os.path.dirname(__file__) + '/progress-circle/progress-circle.js'
        js = open(fn2).read()
        display(HTML(style + circle_html % (js, circle, circle)))
    
    @classmethod
    def make_bar(cls, idx, max_frames, circle):
        '''work with jupyter notebook.
        '''
        percent = str(100*idx/max_frames)
        display(Javascript(js_circle % (circle, percent)))
    
    @classmethod
    def log_progress(cls, sequence, every=1, size=None):
        index = 0
        circle = str(uuid.uuid4())
        cls.init_display(circle)
    
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                cls.make_bar(index, size, circle)
            yield record

class BarProgress(object):
    @classmethod
    def init_display(cls, divid0, divid1, color='#0080FF'):
        my_html = pb.format(id0=divid0, id1=divid1, color=color, msg='0')
        display(HTML(my_html))
    
    @classmethod
    def make_bar(cls, idx, max_frames, divid0, divid1):
        '''work with jupyter notebook.
        '''
        percent = str(100*idx/max_frames)
        s0 = "$('div#{id0}').text({percent});".format(id0=divid0, percent=percent)
        s1 = "$('div#{id1}').width('{percent}%');".format(id1=divid1, percent=percent)
        s = s0 + '\n' + s1
        display(Javascript(s))
    
    @classmethod
    def log_progress(cls, sequence, every=1, size=None, color='#0080FF'):
        index = 0
        divid0 = str(uuid.uuid4())
        divid1 = str(uuid.uuid4())
        cls.init_display(divid0, divid1, color=color)
    
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                cls.make_bar(index, size, divid0, divid1)
            yield record


class ProgressBarTrajectory(object):
    """Simple progress bar/circle for Jupyter notebook

    Parameters
    ----------
    traj : Trajectory-like
    style: str {'bar', 'circle', 'tqdm', callable}, 'bar'
    kwargs : additional keyword arguments

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> from pytraj.utils.progress import ProgressBarTrajectory
    >>> p = ProgressBarTrajectory(traj, style='bar', every=20)
    >>> pt.molsurf(p) # make sure to use Jupyter notebook
    >>> p = ProgressBarTrajectory(traj, style='circle', every=20)
    >>> pt.molsurf(p) # make sure to use Jupyter notebook
    >>> p = ProgressBarTrajectory(traj, style='tqdm')
    >>> pt.molsurf(p) # make sure to use Jupyter notebook
    """
    def __init__(self, traj, style='bar', **kwargs):
        self.traj = traj
        self.style = style
        self.params = kwargs

        for att in dir(traj):
            if not (att.startswith('__') or att == 'xyz'):
                setattr(self, att, getattr(traj, att))
            if att in ['__getstate__', '__setstate__', '_split_iterators']:
                setattr(self, att, getattr(traj, att))

    @property
    def xyz(self):
        # set xyz here to avoid eager evaluation for TrajectoryIterator
        return self.traj.xyz

    def __getitem__(self, index):
        return self.traj[index]

    def __iter__(self):

        if self.style == 'bar':
            every = self.params.get('every', 1)
            color = self.params.get('color', '#0080FF')
            my_iter = BarProgress.log_progress(self.traj, every=every,
                size=self.n_frames, color=color)
        elif self.style == 'circle':
            every = self.params.get('every', 1)
            my_iter = CircleProgress.log_progress(self.traj, every=every,
                size=self.n_frames)
        elif self.style == 'tqdm':
            from tqdm import tqdm_notebook
            my_iter = tqdm_notebook(self.traj, total=self.n_frames)
        elif callable(self.style):
            my_iter = self.style(self.traj, **self.params)
        else:
            raise ValueError("style must be string or callable object")

        for frame in my_iter:
            yield frame
