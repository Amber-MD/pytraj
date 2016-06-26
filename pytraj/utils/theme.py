"""highly experimental. Use with your own risk

    from pytraj.utils import theme
    theme.oceans16() # oceans16 https://github.com/dunovank/jupyter-themes

Retart your notebook to reset to default Jupyter theme

If you want to set global theme for your notebook, it's better to install jupyter-themes
https://github.com/dunovank/jupyter-themes

"""
import os
from IPython import display

style = """
<style>
{}
</style>
"""

def _get_theme(css_file):
    dirname = os.path.dirname(os.path.abspath(__file__))
    css_file = os.path.join(dirname, 'css', css_file)
    css = open(css_file).read()
    return display.HTML(style.format(css))

def oceans16():
    return _get_theme('oceans16.css')
